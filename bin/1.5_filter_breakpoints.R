#!/usr/bin/env Rscript
library(logging)
basicConfig()
library(argparse)

parser <- ArgumentParser()
parser$add_argument("dataset", help = "Path to processed logR data")
parser$add_argument("samplename", help = "Sample to process")
parser$add_argument("segmentation",
                    help = "Merged segmentation file [Delimited format]")
parser$add_argument("metadata",
                    help = paste("Metadata table giving",
                                 "purity/tetraploidy/host details",
                                 "for samples [Delimited format]"))
parser$add_argument("savepath",
                    help = paste("Directory to write results into.",
                                 "Will be created if it does not exist.",
                                 "Contents may be overwritten"))
parser$add_argument("-n", "--segmentation_min_size",
                    type = "integer", default = 10, required = FALSE)
parser$add_argument("-g", "--segmentation_penalty",
                    type = "double", default = 200, required = FALSE)
parser$add_argument("-l", "--filter_tolerance",
                    type = "double", default = 0.25, required = FALSE)

args <- parser$parse_args()

for (argname in names(args)) {
    logdebug("%s = %s", argname, args[[argname]])
}


if (!dir.exists(args$dataset)) {
    stop(sprintf("Dataset %s not found", args$dataset))
}
if (!file.exists(args$segmentation)) {
    stop(sprintf("Could not find segmentation file %s", args$segmentation))
}
if (!dir.exists(args$savepath)) {
    dir.create(args$savepath, recursive = TRUE)
}
if (!dir.exists(args$savepath)) {
    stop("Could not create output path %s", args$savepath)
}
if (!file.exists(args$metadata)) {
    stop("Could not find file %s", args$metadata)
}

suppressPackageStartupMessages({
    library(data.table)
    library(progress)
    library(cnpipe)
    library(arrow)
    library(dplyr)
    library(segmentation)
})


add_segmentation_to_table <- function(data_table, seg) { # nolint start
    chroms <- data_table[, sort(unique(CHROM))]
    data_table[, is_bk := 0]
    data_table[, segment_id := 0]

    for (chrom in chroms) {
        segc <- seg[[chrom]]
        tmp <- data_table[CHROM == chrom & excluded_from_segmentation == FALSE]
        tmp[, I := .I]
        bk_starts <- tmp[I %in% segc$sta, .(CHROM = chrom, START)]
        data_table[bk_starts, is_bk := 1, on = c("CHROM", "START")]
    }

    data_table[, segment_id := cumsum(is_bk), by = CHROM]
    data_table[segment_id == 0, segment_id := 1]
    data_table[, is_bk := NULL] # nolint end
}


add_segment_medians_to_table <- function(data_table) {
    medians <- data_table[excluded_from_segmentation == FALSE,
                          .(segment_median = round(median(total_cn))),
                          by = .(CHROM, segment_id)]
    data_table[medians,
               segment_median := i.segment_median,
               on = c("CHROM", "segment_id")]
}


tost_prefilter <- function(dta, seg, tol = 0.3, alpha = 0.05) { #nolint start
    # Returns TRUE if the distributions a and b are equivalent,
    # within tolerance `tol`
    tost_equivalent <- function(a, b, tol = 0.5, alpha = 0.05) {
        test1 <- t.test(a, b, mu = tol, var.equal = TRUE,
                        alternative = "less")$p.value
        test2 <- t.test(a, b, mu = -tol, var.equal = TRUE,
                        alternative = "greater")$p.value
        test1 < alpha & test2 < alpha
    }

    dta[, bk := rep(seq(1, seg$nIntervals), times = seg$Lengde)]
    dta[, segment_median := median(total_cn), by = bk]

    nseg <- dta[, max(bk)]

    if (nseg > 1) {
        tost_results <- sapply(1:(nseg - 1), function(i) {
            result <- tost_equivalent(dta[segment_id == i, total_cn],
                                      dta[segment_id == (i + 1), total_cn],
                                      tol, alpha)
            return(result)
        })
        approved_breaks <- (2:nseg)[!tost_results]
    } else {
        approved_breaks <- 1
    }

    is_breakpoint_start <- dta[bk %in% approved_breaks,
                               .(START = min(START)),
                               by = bk]
    dta[, tmp := 0]
    dta[is_breakpoint_start, tmp := 1, on = c("START")]
    dta[, filtered_segment_id := cumsum(tmp) + 1]
    dta[, tmpI := .I]

    result <- dta[, .(Lengde = .N, sta = min(tmpI), mean = mean(total_cn)),
                  by = filtered_segment_id]
    result <- list(Lengde = result[, Lengde],
                   sta = result[, sta],
                   mean = result[, mean],
                   nIntervals = dta[, max(filtered_segment_id)])

    dta[, tmp := NULL]
    dta[, tmpI := NULL]
    dta[, unfiltered_segment_id := bk]
    dta[, bk := NULL]
    dta[, segment_median := NULL]

    return(result)
} # nolint end

# Merge adjacent segments where they have equal copy number
merge_equal_copynumber_segments <- function(data_table) { # nolint start
    segments <- data_table[, .(START = min(START)),
                           by = .(segment_id, segment_median, CHROM)]
    segments[, prev_segment_median := shift(segment_median, type = "lag"),
             by = CHROM]
    segments[, is_bk := TRUE]
    selecter <- segments[segment_median == prev_segment_median, is_bk := FALSE]
    data_table[, is_bk := 0]
    data_table[, segment_id := 0]
    data_table[selecter[is_bk == TRUE], is_bk := 1, on = c("CHROM", "START")]
    data_table[, segment_id := cumsum(is_bk), by = CHROM]
    data_table[segment_id == 0, segment_id := 1]
    data_table[, is_bk := NULL]
} # nolint end


loginfo("Loading data")
metadata <- fread(args$metadata, na.strings = "")
merged_breakpoints <- fread(args$segmentation)
ds <- arrow::open_dataset(args$dataset)
ds %>% filter(samplename == args$samplename) %>% collect() -> dt
setDT(dt)
dt[, CHROM := as.character(CHROM)]
setorder(dt, CHROM, START, END, samplename)

loginfo("Segmenting")
chroms <- dt[, sort(unique(CHROM))]
seg <- lapply(chroms, function(chrom) {
    available_breakpoints <- merged_breakpoints[CHROM == chrom, start]
    cn_data <- dt[excluded_from_segmentation == FALSE & CHROM == chrom,
                  total_cn]
    restricted_pcf(cn_data, available_breakpoints,
                   args$segmentation_min_size,
                   args$segmentation_penalty)
})
names(seg) <- chroms
invisible({
    add_segmentation_to_table(dt, seg)
    add_segment_medians_to_table(dt)
})

loginfo("Filtering segmentation")
segf <- lapply(chroms, function(chrom) {
    print(chrom)
    suppressWarnings(tost_prefilter(dt[(excluded_from_segmentation == FALSE &
                                        CHROM == chrom)],
                                    seg[[chrom]],
                                    tol = args$filter_tolerance,
                                    alpha = 0.05))
})
names(segf) <- chroms
invisible({
    add_segmentation_to_table(dt, segf)
    add_segment_medians_to_table(dt)
})

loginfo("Merging equal copy number segments")
invisible({
    merge_equal_copynumber_segments(dt)
    add_segment_medians_to_table(dt)
})

outdir <- file.path(args$savepath, "dataset")
loginfo("Writing output to %s", outdir)
setorder(dt, CHROM, START, END, samplename)
arrow::write_dataset(dt, outdir, partitioning = "samplename")

clone_ <- metadata[tumour == dt[1, samplename], clone]
all_breakpoints <- dt[, .(START = min(START), is_bk = TRUE),
                      by = .(CHROM, segment_id)]
all_breakpoints[, clone := clone_]
outfile <- file.path(args$savepath,
                     paste0(dt[1, samplename], ".all_breakpoints.csv"))
loginfo("Writing breakpoint summary to %s", outfile)
fwrite(all_breakpoints, outfile)
