#!/usr/bin/env Rscript

library(argparse)
library(data.table)

parser <- ArgumentParser()
parser$add_argument("--input_path",
    help = "Path to the directory containing the cleaned read count files")
parser$add_argument("--output_file",
    help = "Path to the output file")
args <- parser$parse_args()

if(!dir.exists(args$input_path)) {
    stop("Input path does not exist")
}

files <- list.files(path = args$input_path,
    pattern = "*.cleaned_counts.tsv.gz",
    full.names = TRUE)

read_count_file <- function(filename) {
    dt <- fread(filename)
    samplename <- dt[1, samplename]
    setnames(dt, old = "COUNT", new = samplename)
    setkey(dt, CHROM, START, END)
    dt[, 1:4]
}

merge_count_file <- function(dt, filename) {
    dt2 <- read_count_file(filename)
    dt <- merge(dt, dt2, all = TRUE)
    dt
}

merged <- Reduce(merge_count_file,
    files[2:length(files)],
    init = read_count_file(files[1]))

fwrite(merged, args$output_file, sep = "\t", na = "N/A")
