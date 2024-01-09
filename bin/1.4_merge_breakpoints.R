#!/usr/bin/env Rscript
library(logging); basicConfig()
library(argparse)

parser = ArgumentParser()
parser$add_argument("inputpath",
                    help = "Path to segmentation files [Delimited format]")
parser$add_argument("metadata",
                    help = paste("Metadata table giving purity",
                                 "/tetraploidy/host details for samples ",
                                 "[Delimited format]", sep = ""))
parser$add_argument("clone",
                    help = "The clone for which to merge segments")
parser$add_argument("savepath",
                    help = paste("Directory to write results into. Will be",
                                 "created if it does not exist. Contents may",
                                 "be overwritten"))
parser$add_argument("-k", "--merge_window_size", type = "integer",
                    help = "Minimum size of a segment after merging [min=1]",
                    default = 10)

args <- parser$parse_args()

for (argname in names(args)) {
    loginfo("%s = %s", argname, args[[argname]])
}


if (!dir.exists(args$inputpath)) {
    stop(sprintf("Could not find input path %s", args$inputpath))
}
if (!dir.exists(args$savepath)) {
    dir.create(args$savepath, recursive=TRUE)
}
if (!dir.exists(args$savepath)) {
    stop("Could not create output path %s", args$savepath)
}
if (!file.exists(args$metadata)) {
    stop("Could not find file %s", args$metadata)
}
if (args$merge_window_size < 1) {
    stop("Merge window size is below minimum allowed value (1)")
}

library(data.table)
metadata <- fread(args$metadata, na.strings="")
clone_samples <- metadata[!is.na(tumour) & clone == args$clone, tumour]
segmentation_files <- file.path(args$inputpath,
                                paste(clone_samples,
                                      "filtered_segmentation",
                                      "csv",
                                      sep = "."))
if (!all(file.exists(segmentation_files))) {
    stop(sprintf("Missing input file - %s",
                 segmentation_files[!file.exists(segmentation_files)][1]))
}

get_merged_breakpoints <- function(segmentations_dt, cutoff,
                                   use_weighted_mean = TRUE) {

    assign_merged_segment_id <- function(starts, cutoff) {
        if (length(starts) == 1) {
            return (as.integer(1))
        } else {
            as.integer(cutree(hclust(dist(starts)^2, method = "centroid"),
                              h = cutoff^2))
        }
    }

    chrom_maxes <- segmentations_dt[, .(end=max(end)), by = CHROM]
    unique_starts <- unique(segmentations_dt, by = c("CHROM", "start"))
    setorder(unique_starts, CHROM, start)
    merged <- unique_starts[,
                            .(start,
                              segmentID = assign_merged_segment_id(start,
                                                                   cutoff)),
                            by = CHROM]
    if (use_weighted_mean) {
        tmp <- copy(segmentations_dt)
        merged <- tmp[merged, segmentID := i.segmentID,
                      on = c("CHROM", "start")]
    }
    merged <- merged[, .(start=floor(mean(start))),
                     by = .(CHROM, segmentID)][order(CHROM, segmentID)]

    # Averaging may cause some segments to become smaller than the cutoff
    niter <- 0
    while (TRUE) {
        niter <- niter + 1
        merged[, diff := start - shift(start), by = CHROM]
        merged[diff < cutoff, start := start + cutoff - diff]
        success <- merged[diff < cutoff, .N] == 0
        merged[, diff := NULL]
        if (success) {
            break
        }
        if (niter == 10) {
            logwarn(paste("Some segments remain smaller than the",
                    "cutoff after %d attempts to fix this"),
            niter)
            break
        }
    }

    merged[, end := data.table::shift(start - 1, type="lead"), by = CHROM]

    # Set the end breakpoint to the end of each chromosome
    merged[merged[, .(start=max(start)), by = CHROM][chrom_maxes, ,
                                                     on = "CHROM"],
           end := i.end,
           on = c("start", "CHROM")]

    merged[, K := end - start + 1]
    return (merged)
}



segmentation <- rbindlist(lapply(segmentation_files, fread))
merged_breakpoints <- get_merged_breakpoints(copy(segmentation),
                                             args$merge_window_size,
                                             use_weighted_mean = TRUE)

output_file <- file.path(args$savepath,
                         sprintf("merged_breakpoints_%s.csv", args$clone))
fwrite(merged_breakpoints, output_file)

loginfo("Output written to %s", output_file)
