#!/usr/bin/env Rscript
library(data.table)
library(magrittr)
library(argparse)
parser <- ArgumentParser(
    description = "Combine ploidy estimates from multiple samples")
parser$add_argument("-i", "--input",
    help = "Input files to combine",
    nargs = "+",
    required = TRUE)
parser$add_argument("-o", "--output",
    help = "Output file",
    required = TRUE)

args <- parser$parse_args()

lapply(args$input, function(x) {
    if (file.exists(x)) {
        fread(x)
    } else {
        stop(paste0("File ", x, " does not exist"))
    }
}) %>%
  rbindlist() %>%
  fwrite(args$output)
