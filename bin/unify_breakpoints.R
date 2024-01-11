#!/usr/bin/env Rscript
library(logging)

suppressPackageStartupMessages({
    library(data.table)
})

segmentation_files <- list.files(".", "*all_breakpoints.csv")

all_breakpoints <- rbindlist(lapply(file.path(".",
                                              segmentation_files), fread))
unified_breakpoints <- unique(all_breakpoints,
                              by = c("CHROM", "START", "clone"))[order(clone,
                                                                       CHROM,
                                                                       START)]

fwrite(unified_breakpoints, "unified_breakpoints.csv")
