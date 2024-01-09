#!/usr/bin/env Rscript
# 1 - Calculate the total copy number for each sample
# 2 - Segment the whole genome
# 3 - Filter the segments

# Handle command line arguments before loading any libraries
library(logging)
basicConfig()

library(argparse)
parser <- ArgumentParser()

parser$add_argument("dataset", help = "Path to processed logR data")
parser$add_argument("samplename", help = "Sample to process")
parser$add_argument("metadata",
                    help = paste("Metadata table giving",
                                 "purity/tetraploidy/host details for samples",
                                 "[Delimited format]"))
parser$add_argument("savepath",
                    help = paste("Directory to write results into. Will be",
                                 "created if it does not exist. Contents may",
                                 "be overwritten"))
parser$add_argument("-e", "--excluded_bins",
                    help = paste("Table of bins to exclude from segmentation",
                                 "[Delimited format]"))
parser$add_argument("-d", "--diploid_ploidy_value",
                    type = "double", default = 2.0, required = FALSE)
parser$add_argument("-t", "--tetraploid_ploidy_value",
                    type = "double", default = 4.0, required = FALSE)
parser$add_argument("-p", "--pseudoautosome_end",
                    type = "integer", default = 0, required = FALSE)
parser$add_argument("-n", "--segmentation_min_size",
                    type = "integer", default = 10, required = FALSE)
parser$add_argument("-g", "--segmentation_penalty",
                    type = "double", default = 200, required = FALSE)
parser$add_argument("-w", "--winsorize", action = "store_true",
                    help = "Winsorize copy number values before segmentation")
parser$add_argument("-k", "--winsorization_window_size",
                    type = "double", default = 25, required = FALSE)
parser$add_argument("-y", "--winsorization_strength",
                    type = "double", default = 2.5, required = FALSE)
parser$add_argument("-l", "--filter_tolerance",
                    type = "double", default = 0.25, required = FALSE)
parser$add_argument("-a", "--do_ploidy_assessment",
                    action = "store_true", default = FALSE, required = FALSE)
parser$add_argument("-c", "--consolidate_logr",
                    action = "store_true", default = FALSE, required = FALSE,
                    help = paste("Select this flag if logr is divided into",
                                 "separate columns for",
                                 "T_logr_denoised_male/female"))

args <- parser$parse_args()

for (argname in names(args)) {
    loginfo("%s = %s", argname, args[[argname]])
}

if (!dir.exists(args$savepath)) dir.create(args$savepath, recursive = TRUE)
if (!dir.exists(args$savepath)) stop("Unable to create output directory")

if (!dir.exists(args$dataset)) {
    stop(sprintf("Dataset %s not found", args$dataset))
}

if (!file.exists(args$metadata)) {
    stop(sprintf("File %s not found", args$metadata))
}

if (!is.null(args$excluded_bins)) {
    if (!file.exists(args$excluded_bins)) {
        stop(sprintf("File %s not found", args$excluded_bins))
    }
}

suppressPackageStartupMessages({
    library(cnpipe)
    library(arrow)
    library(dplyr)
    library(data.table)
    library(segmentation)
})


load_table_from_dataset <- function(path, samplename, bins_to_exclude,
                                    metadata, use_standardized = FALSE,
                                    pseudoautosome_end = 0) { # nolint start
    ds <- arrow::open_dataset(path)
    samplename_ <- samplename
    ds %>% filter(samplename == samplename_) %>% collect() -> data_table
    setDT(data_table)
    data_table[, CHROM := as.character(CHROM)]
    data_table[, START := as.integer(round(START, -3))]
    setorder(data_table, CHROM, START, END, samplename)
    data_table[, Index := seq(1, .N), by = .(CHROM)]
    data_table[, excluded_from_segmentation := FALSE]
    data_table[bins_to_exclude, excluded_from_segmentation := TRUE,
               on = c("CHROM", "START", "END")]
    data_table[START %% 1000 == 0, START := START+1]
    if (args$consolidate_logr) {
        data_table <- consolidate_logr(data_table, metadata, use_standardized)
    }
    data_table <- add_expected_host_copy_number_to_table(data_table,
                                                         metadata,
                                                         pseudoautosome_end)
    setcolorder(data_table, c("CHROM", "START", "END", "samplename"))
    return(data_table)
} # nolint end


#' Constructs a single reference tumour logR column, replacing the multiple
#' standardized/denoised, male/female/neutral columns that exist in the input
consolidate_logr <- function(data_table, metadata, use_standardized=FALSE) { # nolint start
    tumour_ <- data_table[1, samplename]
    if (startsWith(tumour_, "s")) {
        logdebug("Samplename starts with 's', which may be an older pipeline artifact")
    }
    host_sex <- metadata[tumour == tumour_, hostSex]
    if (use_standardized) {
        data_table[CHROM != "X" & CHROM != "Y", T_logr := T_logr_standardized]
        if (host_sex == "M" | host_sex == "N") {
            data_table[CHROM == "X" | CHROM == "Y",
                       T_logr := T_logr_standardized_male]
        } else {
            data_table[CHROM == "X" | CHROM == "Y",
                       T_logr := T_logr_standardized_female]
        }
    } else {
        data_table[CHROM != "X" & CHROM != "Y", T_logr := T_logr_denoised]
        if (host_sex == "M" | host_sex == "N") {
            data_table[CHROM == "X" | CHROM == "Y",
                       T_logr := T_logr_denoised_male]
        } else {
            data_table[CHROM == "X" | CHROM == "Y",
                       T_logr := T_logr_denoised_female]
        }
    }
    data_table[, c("T_logr_denoised",
                   "T_logr_denoised_male",
                   "T_logr_denoised_female",
                   "T_logr_standardized",
                   "T_logr_standardized_male",
                   "T_logr_standardized_female") := NULL]
    return(data_table)
} # nolint end


add_expected_host_copy_number_to_table <- function(data_table,
                                                   metadata,
                                                   pseudoautosome_end=0) { # nolint start
    tumour_ <- data_table[1, samplename]
    if (startsWith(tumour_, "s")) {
        logdebug("Samplename starts with 's', which may be an older pipeline artifact")
    }
    host_sex <- metadata[tumour == tumour_, hostSex]

    data_table[, expected_host_cn := NA_integer_]
    data_table[CHROM != "X" & CHROM != "Y",
               expected_host_cn := 2]
    data_table[CHROM == "X" & END <= pseudoautosome_end,
               expected_host_cn := 2]
    data_table[CHROM == "X" & START > pseudoautosome_end,
               expected_host_cn := ifelse(host_sex == "F", 2, 1)]
    data_table[CHROM == "Y",
               expected_host_cn := ifelse(host_sex == "F", 0, 1)]
} # nolint end


add_copy_number_to_table <- function(data_table, metadata,
                                     diploid_ploidy_value,
                                     tetraploid_ploidy_value,
                                     winsorize = FALSE,
                                     tau = 2.5, k = 25) { # nolint start
    stopifnot("expected_host_cn" %in% colnames(data_table))
    stopifnot("T_logr" %in% colnames(data_table))

    tumour_ <- data_table[1, samplename]
    if (startsWith(tumour_, "s")) {
        logdebug("Samplename starts with 's', which may be an older pipeline artifact")
    }
    purity <- metadata[tumour == tumour_, purity]
    tetraploid <- metadata[tumour == tumour_, isTetraploid]
    host_sex <- metadata[tumour == tumour_, hostSex]

    # For dogs, the pseudoautosome region on X ends around 6.645e6
    data_table[, total_cn := 0]

    if (tetraploid) {
        data_table[, total_cn := calc_total_copynumber(T_logr, purity,
                                                       tetraploid_ploidy_value,
                                                       expected_host_cn, 2.0)]
        if (winsorize) {
            data_table[, total_cn := winsorize(total_cn, tau, k)]
        }
    } else {
        data_table[, total_cn := calc_total_copynumber(T_logr, purity,
                                                       diploid_ploidy_value,
                                                       expected_host_cn, 2.0)]
        if (winsorize) {
            data_table[, total_cn := winsorize(total_cn, tau, k)]
        }
    }
} # nolint end


make_segmentation <- function(data_table, kmin=10, gamma=100) { # nolint start
    chroms <- data_table[, sort(unique(CHROM))]

    segmentations <- lapply(chroms, function(chrom) {
                                print(chrom)
        pcf(data_table[CHROM == chrom, total_cn], kmin, gamma)
    })

    names(segmentations) <- chroms
    return(segmentations)
} # nolint end


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
    data_table[, is_bk := NULL]
} # nolint end


add_segment_medians_to_table <- function(data_table) {
    medians <- data_table[excluded_from_segmentation == FALSE,
                          .(segment_median = round(median(total_cn))),
                          by = .(CHROM, segment_id)]
    data_table[medians, segment_median := i.segment_median,
               on = c("CHROM", "segment_id")]
}


assess_ploidy <- function(data_table) {
    data_table[excluded_from_segmentation == FALSE,
               sum(sign(total_cn - segment_median))]
}


assess_ploidy_range <- function(data_table, metadata) { # nolint start
    ploidies <- seq(3.44, 4.26, 0.01)
    results <- sapply(ploidies, function(ploidy) {
        add_copy_number_to_table(data_table, metadata, ploidy / 2, ploidy)
        add_segment_medians_to_table(data_table)
        result <- assess_ploidy(data_table)
        print(sprintf("%.2f %.3f", ploidy, data_table[, mean(segment_median)]))
        return(result)
    })
    return(data.table(diploidy = ploidies / 2, tetraploidy = ploidies,
                      results))
} # nolint end


tost_prefilter <- function(dta, seg, tol=0.3, alpha=0.05) { # nolint start
    # Returns TRUE if the distributions a and b are equivalent within
    # tolerance `tol`
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
                              .(START = min(START)), by = bk]
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


convert_segmentation_to_data_table <- function(segm, samplename_) {
    seg_dt <- rbindlist(lapply(chroms, function(chrom) {
        segmc <- segm[[chrom]]
        data.table(CHROM = chrom,
                   start = segmc$sta,
                   end = segmc$sta + segmc$Lengde - 1,
                   K = segmc$Lengde,
                   mean = segmc$mean)
    }))
    seg_dt[, samplename := samplename_]
    return(seg_dt)
}


check_datatable_has_columns <- function(data_table, required_columns) {
    if (!all(required_columns %in% colnames(data_table))) {
        stop(sprintf("Metadata table is missing at least one required column (%s)",
                     paste(required_columns, collapse = ",")))
    }
}

metadata <- fread(args$metadata, na.strings = "")
check_datatable_has_columns(metadata, c("tumour", "hostSex", "purity",
                                        "isTetraploid"))

bins_to_exclude <- fread(args$excluded_bins)
if (is.null(bins_to_exclude)) {
    bins_to_exclude <- data.table(CHROM = vector("character"),
                                  START = vector("integer"),
                                  END = vector("integer"))
}
bins_to_exclude[, CHROM := as.character(CHROM)]
bins_to_exclude[, START := as.integer(round(START, -3))]

loginfo("Loading %s", args$filename)
dt <- load_table_from_dataset(args$dataset, args$samplename, bins_to_exclude,
                              metadata, use_standardized = FALSE,
                              pseudoautosome_end = args$pseudoautosome_end)
check_datatable_has_columns(dt, "T_logr")
samplename_ <- dt[1, samplename]
add_copy_number_to_table(dt, metadata,
                         args$diploid_ploidy_value,
                         args$tetraploid_ploidy_value,
                         args$winsorize,
                         args$winsorization_strength,
                         args$winsorization_window_size)

loginfo("Segmenting")
seg <- make_segmentation(dt[excluded_from_segmentation == FALSE],
                         args$segmentation_min_size,
                         args$segmentation_penalty)
add_segmentation_to_table(dt, seg)
add_segment_medians_to_table(dt)
loginfo("(%s) (%.2f,%.2f) TOTAL_SEGMENTS: %d, PLOIDY SCORE: %d",
        samplename_,
        args$diploid_ploidy_value,
        args$tetraploid_ploidy_value,
        sum(sapply(seg, `[[`, "nIntervals")),
        assess_ploidy(dt))

chroms <- dt[, sort(unique(CHROM))]
loginfo("Filtering segmentation")
segf <- lapply(chroms, function(chrom) {
    suppressWarnings(tost_prefilter(dt[(excluded_from_segmentation == FALSE &
                                        CHROM == chrom)],
                                    seg[[chrom]], tol = args$filter_tolerance,
                                    alpha = 0.05))
})
names(segf) <- chroms

add_segmentation_to_table(dt, segf)
add_segment_medians_to_table(dt)
loginfo("(%s) (%.2f,%.2f) TOTAL_SEGMENTS: %d, PLOIDY SCORE: %d",
        samplename_,
        args$diploid_ploidy_value,
        args$tetraploid_ploidy_value,
        sum(sapply(segf, `[[`, "nIntervals")),
        assess_ploidy(dt))

segmentation <- convert_segmentation_to_data_table(seg, samplename_)
filtered_segmentation <- convert_segmentation_to_data_table(segf, samplename_)

outpath <- file.path(args$savepath, "dataset")
loginfo("Writing data to %s", outpath)
setorder(dt, CHROM, START, END, samplename)
arrow::write_dataset(dt, outpath, partitioning = "samplename")

fwrite(segmentation,
       file.path(args$savepath, paste0(samplename_, ".segmentation.csv")))
fwrite(filtered_segmentation,
       file.path(args$savepath,
                 paste0(samplename_, ".filtered_segmentation.csv")))

if (args$do_ploidy_assessment) {
    loginfo("Running ploidy assessment")
    ploidy_assessment <- assess_ploidy_range(dt, metadata)
    ploidy_assessment_filename <- file.path(args$savepath,
                                            paste0(samplename_,
                                                   ".ploidy_assessment.csv"))
    loginfo("Writing ploidy assessment result to %s",
            ploidy_assessment_filename)
    fwrite(ploidy_assessment, ploidy_assessment_filename)

    # Estimate the best-fit diploidy
    minimum_sq_score <- ploidy_assessment[, (which.min((results / 1e5) ^ 2))]
    range_to_interpolate <- seq(max(1, minimum_sq_score - 3),
                                min(nrow(ploidy_assessment),
                                    minimum_sq_score + 3))
    linmod <- lm(results ~ diploidy,
                 data = ploidy_assessment[range_to_interpolate])
    diploidy_estimate <- -coef(linmod)["(Intercept)"] / coef(linmod)["diploidy"]
    tetraploidy_estimate <- 2 * diploidy_estimate
    ploidy_estimate <- if (metadata[tumour == samplename_, isTetraploid]) {
        tetraploidy_estimate
    } else {
        diploidy_estimate
    }
    loginfo("Estimated ploidy: %.3f", ploidy_estimate)

    # Write some updated output files
    new_row <- metadata[tumour == samplename_]
    new_row[, estimatedPloidy := ploidy_estimate]
    fwrite(new_row[tumour == samplename_],
           file.path(args$savepath, paste0(samplename_, ".metadata.csv")))

    add_copy_number_to_table(dt, metadata, diploidy_estimate,
                             tetraploidy_estimate,
                             args$winsorize,
                             args$winsorization_strength,
                             args$winsorization_window_size)
    add_segment_medians_to_table(dt)
}
