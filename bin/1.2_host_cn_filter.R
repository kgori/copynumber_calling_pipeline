#!/usr/bin/env Rscript
library(argparse)
library(logging); basicConfig()
parser <- ArgumentParser()
parser$add_argument("-i", "--inputpath", required = TRUE,
                    help = "Path to directory of input files")
parser$add_argument("-o", "--outputpath", required = TRUE,
                    help = "Path to directory into which to write output")
parser$add_argument("-m", "--metadata", required = TRUE,
                    help = "Path to metadata file")
parser$add_argument("--modelfile", required = TRUE,
                    help = "Filtering model specification file")
parser$add_argument("--modeltype", default = "RF", choices = c("RF", "NN"),
                    help = paste("Which model to use for prediction.",
                                 "The choices are RF (random forest),",
                                 "or NN (neural network). Default is",
                                 "RF"))
parser$add_argument("--cutoff", default = 0.85, type = "double",
                    help = paste("Probability cutoff for the filter.",
                                 "The filter assesses the probability",
                                 "that a segment is normal diploid in all",
                                 "samples. Segments with a probability",
                                 "less than the cutoff will be ignored",
                                 "during copy number segmentation. Default",
                                 "is 0.85"))
args <- parser$parse_args()

for (argname in names(args)) {
    loginfo("%s = %s (%s)", argname, args[[argname]], class(args[[argname]]))
}

if (!file.exists(args$metadata)) {
    stop(sprintf("Metadata table %s could not be found", args$metadata))
}

if (!dir.exists(args$inputpath)) {
    stop(sprintf("Input directory %s could not be found", args$inputpath))
}

if (!file.exists(args$modelfile)) {
    stop(sprintf("Model file %s could not be found", args$modelfile))
}

if (!grepl(args$modeltype, args$modelfile)) {
    logwarn(sprintf("Model file %s may not match model type %s",
                 args$modelfile, args$modeltype))
}

if (!dir.exists(args$outputpath)) {
    dir.create(args$outputpath, recursive = TRUE)
}

if (!dir.exists(args$outputpath)) {
    stop(sprintf("Output directory %s could not be found or created",
                 args$outputpath))
}

suppressPackageStartupMessages({
library(data.table)
library(arrow)
library(dplyr)
library(randomForest)
library(segmentation)
library(cnpipe)
library(naturalsort)
})

# Get current path of this executing script
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name,
                                                           initial.options)])
script.basename <- normalizePath(dirname(script.name))
if (identical(script.basename, character(0))) {
    script.basename <- normalizePath(".")
}
logdebug("Script execution directory = %s", script.basename)

filtermodel <- readRDS(args$modelfile)

# Alias the model to the appropriate name
if (args$modeltype == "RF") {
    forest <- filtermodel
} else {
    nn <- filtermodel
}

#' Takes a vector of input data, and returns a summary of its distribution
#' as a vector of length 128 (which matches the neural network)
to_nn_input <- function(datavec) {
    y <- density(datavec, n = 128, from = 0, to = 5)$y
    y / max(y)
}

#' Make a prediction of whether the data passes the filter,
#' using the trained neural network
nn_predict <- function(data, network) {
    suppressMessages(attach(network))
    on.exit(detach(network))
    stopifnot(length(data) == nrow(network$w1))
    relu <- network$relu
    network$predict(data)[1, 1]
}


metadata <- fread(args$metadata, na.strings = c("NA", "na", "N/A", "n/a", ""))
ds <- arrow::open_dataset(file.path(args$inputpath, "dataset"))


hosts_to_exclude_from_panel <-
    if ("excludedFromPanel" %in% colnames(metadata)) {
        metadata[excludedFromPanel == TRUE, host]
    } else {
        c()
    }

all_hosts <- setdiff(metadata[!is.na(tumour) & !is.na(host), host],
                     hosts_to_exclude_from_panel)
female_hosts <- setdiff(metadata[!is.na(tumour) &
                        !is.na(host) &
                        hostSex == "F", host],
                        hosts_to_exclude_from_panel)

# filenames <- list.files(args$inputpath, pattern = "*.fst",
#                         include.dirs = TRUE, full.names = TRUE)
# file.exists(filenames)

# Build a matrix of host copy number estimates
tumours_with_matched_hosts <- unique(metadata[data.table(host = all_hosts), ,
                                     on = "host"], by = "host")[, tumour]
# tumours_with_matched_hosts_filenames <-
#     sapply(tumours_with_matched_hosts,
#         function(tumour) {
#             grep(tumour, filenames, value = T)
#         })

loginfo("Loading host copy number estimates")
ds %>%
    filter(samplename %in% tumours_with_matched_hosts) %>%
    select(c("CHROM", "START", "END", "hostname", "host_cn")) %>%
    collect() -> dt

dt <- dcast(dt, CHROM + START + END ~ hostname, value.var = "host_cn")
invisible(dt[, lapply(.SD, winsorize), .SDcols = all_hosts])


loginfo("Segmenting host copy numbers")

# Preselect some locations as likely breakpoints, using convolution
# with a sawtooth kernel
mark_for_segmentation <- function(m) {
    f1 <- c(1:15 / 15, rev(-(1:15 / 15)))
    f2 <- c(1:3 / 3, rev(-(1:3 / 3)))
    c1 <- rowSums(apply(m, 2,
                        function(col) abs(segmentation:::convolve_(col, f1))))
    c2 <- rowSums(apply(m, 2,
                        function(col) abs(segmentation:::convolve_(col, f2))))
    l1 <- quantile(c1, 0.75)
    l2 <- quantile(c2, 0.9)

    which((c1 > l1) | (c2 > l2))
}

chroms <- setdiff(dt[, naturalsort(unique(CHROM))], "Y")
seg <- lapply(chroms, function(chrom) {
    loginfo(" -- Segmenting %s", chrom)
    if (chrom == "X") {
        data <- as.matrix(dt[CHROM == chrom, .SD, .SDcols = female_hosts])
    } else {
        data <- as.matrix(dt[CHROM == chrom, .SD, .SDcols = all_hosts])
    }
    mark <- mark_for_segmentation(data)
    expanding_fast_multipcf(data, mark, kmin = 10, gamma = 7)
})
names(seg) <- chroms

# Cache outputs
saveRDS(dt, file.path(args$outputpath,
                     "host_copy_number_matrix.RDS"))
saveRDS(seg, file.path(args$outputpath,
                       "host_copy_number_segmentation.RDS"))

# Running data through the pretrained NN and RF models

# 1: Convert inputs
segdt <- rbindlist(lapply(names(seg), function(chrom) {
    data.table(CHROM = chrom,
               startpos = seg[[chrom]]$starts,
               endpos = seg[[chrom]]$ends)
}))
segdt[, segID := seq_len(.N), by = CHROM]

dt[, c("index", "startpos", "endpos") := seq_len(.N), by = CHROM]
setkey(dt, CHROM, startpos, endpos)
setkey(segdt, CHROM, startpos, endpos)
ol <- foverlaps(dt, segdt)
setkey(ol, CHROM, index)
setkey(dt, CHROM, index)
dt[ol, segment_index := segID]
dt[, startpos := NULL]
dt[, endpos := NULL]

to_rf_input <- function(m) {
    y <- to_nn_input(m)
    y <- as.data.frame(t(y))
    colnames(y) <- paste("RF", seq_len(ncol(y)), sep = ".")
    y
}

loginfo("Building prediction matrix")
if (args$modeltype == "RF") {
    pred_input <- dt[CHROM != "X" & CHROM != "Y",
                     to_rf_input(as.matrix(.SD)),
                     by = .(CHROM, segment_index),
                     .SDcols = all_hosts]

    if (dt[CHROM == "X", .N] > 0) {
        pred_input <- rbind(pred_input,
                            dt[CHROM == "X",
                               to_rf_input(as.matrix(.SD)),
                               by = .(CHROM, segment_index),
                               .SDcols = female_hosts])
    }
    if (dt[CHROM == "Y", .N] > 0) {
        pred_input <- rbind(pred_input,
                            dt[CHROM == "Y",
                               to_rf_input(as.matrix(.SD)),
                               by = .(CHROM, segment_index),
                               .SDcols = setdiff(all_hosts, female_hosts)])
    }

    old_feature_columns <- setdiff(colnames(pred_input), c("CHROM", "segment_index"))
    new_feature_columns <- attr(forest$terms, "term.labels")
    setnames(pred_input,
             old = old_feature_columns,
             new = new_feature_columns)
} else {
    pred_input <- dt[CHROM != "X" & CHROM != "Y",
                     to_nn_input(as.matrix(.SD)),
                     by = .(CHROM, segment_index),
                     .SDcols = all_hosts]

    if (dt[CHROM == "X", .N] > 0) {
        pred_input <- rbind(pred_input,
                            dt[CHROM == "X",
                               to_nn_input(as.matrix(.SD)),
                               by = .(CHROM, segment_index),
                               .SDcols = female_hosts])
    }

    if (dt[CHROM == "Y", .N] > 0) {
        pred_input <- rbind(pred_input,
                            dt[CHROM == "Y",
                               to_nn_input(as.matrix(.SD)),
                               by = .(CHROM, segment_index),
                               .SDcols = setdiff(all_hosts, female_hosts)])
    }
}

if (args$modeltype == "RF") {
    loginfo("Running random forest filter predictions")
    model_predictions <- predict(forest, pred_input, type = "prob")[, "1"]
} else {
    loginfo("Running neural network filter predictions")
    model_predictions <- pred_input[, .(prediction = nn_predict(V1, nn)), by = .(CHROM, segment_index)]$prediction
}
predictions <- dt[, .(START = min(START), END = max(END)),
                  by = .(CHROM, segment_index)]
predictions[segdt,
            c("START_INDEX", "END_INDEX") := .(i.startpos, i.endpos),
            on = c("CHROM", "segment_index" = "segID")]
predictions[, prediction := model_predictions]


saveRDS(predictions, file.path(args$outputpath,
                               "host_copy_number_predictions.RDS"))


exclusions <- predictions[prediction < args$cutoff,
                               .(CHROM, segment_index)]
bins_to_exclude <- dt[exclusions, .(CHROM, START, END),
                      on = c("CHROM", "segment_index")] 
                  
bins_to_exclude[, START := as.integer(round(START, -3))]
bins_to_exclude[, END := as.integer(END)]

loginfo("Number of bins to exclude: %d", bins_to_exclude[, .N])

outfile <- paste0("host_coverage_", args$modeltype, "_filter_bins_to_exclude.csv")
loginfo("Writing exclusion list to %s", file.path(args$outputpath, outfile))
fwrite(bins_to_exclude,
       file.path(args$outputpath,
                 outfile))
