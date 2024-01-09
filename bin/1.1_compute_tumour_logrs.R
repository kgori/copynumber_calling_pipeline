#!/usr/bin/env Rscript
# We're going to compute the logR for all our Tumours' binned read counts,
# relative to a panel made up of all our Hosts' binned read counts.
# For good measure, and also to deal with the haploidy of the X chromosome,
# we'll compute the logR independently for panels made up of only female and
# only male subsets of the Hosts.

library(argparse)
library(logging)
basicConfig()
parser <- ArgumentParser()
parser$add_argument("-r", "--readcounts", required = TRUE,
                    help = paste("Path to csv file of per-bin readcounts for",
                                 "all tumour and host samples that will be",
                                 "needed"))
parser$add_argument("-o", "--outputpath", required = TRUE,
                    help = paste("Directory into which to write output",
                                 "(will be created if necessary)"))
parser$add_argument("-m", "--metadata", required = TRUE,
                    help = "Path to csv metadata file")
parser$add_argument("--gc-correction", required = FALSE, dest = "gccorrection",
                    default = "", help = "Path to GC interval contents BED file")
parser$add_argument("-n", "--number-of-singular-values", default = 3,
                    type = "integer", dest = "nsvs",
                    help = paste("Number of singular values to use for",
                                 "denoising. Default is 3."))
args <- parser$parse_args()

for (argname in names(args)) {
    loginfo("%s = %s", argname, args[[argname]])
}

if (!file.exists(args$readcounts)) {
    stop(sprintf("Input readcounts file %s not found", args$readcounts))
}

if (!file.exists(args$metadata)) {
    stop(sprintf("Metadata table %s could not be found", args$metadata))
}

if (!dir.exists(args$outputpath)) {
    dir.create(args$outputpath, recursive = TRUE)
}

if (!dir.exists(args$outputpath)) {
    stop(sprintf("Output directory %s could not be found or created",
                 args$outputpath))
}

if (args$gccorrection != "") {
    if (!file.exists(args$gccorrection)) {
        stop(sprintf("GC correction file %s could not be found",
                     args$gccorrection))
    }
}

suppressPackageStartupMessages({
    library(data.table)
    library(cnpipe)
    library(matrixStats)
    library(segmentation)
    library(arrow)
})

loginfo("Loading readcounts from %s", args$readcounts)
readcounts <- fread(args$readcounts)
setkey(readcounts, CHROM, START, END)
readcounts <- readcounts[CHROM %in% c(1:99, "X", "Y")]
logwarn(paste("Ensure that START value is 1-based (i.e. intervals run from",
              "10001-11000, not 10000-11000, for example)"))


loginfo("Loading metadata from %s", args$metadata)
metadata <- fread(args$metadata, na.strings = c("NA", "na", "N/A", "n/a", ""))

if (!("tumour" %in% colnames(metadata))) {
    logerror("Metadata table must contain a column named 'tumour'")
    stop()
}

if (!("host" %in% colnames(metadata))) {
    logerror("Metadata table must contain a column named 'host'")
    stop()
}

if (!("hostSex" %in% colnames(metadata))) {
    logerror("Metadata table must contain a column named 'hostSex'")
    stop()
}

if (!("excludedFromPanel" %in% colnames(metadata))) {
    logerror("Metadata table does not contain a column named 'excludedFromPanel'")
    stop()
}

if (!("tumourContaminated" %in% colnames(metadata))) {
    logerror("Metadata table does not contain a column named 'tumourContaminated'")
    stop()
}

loginfo("Checking for unrecognised samples")
# Do some cross-referencing of samples
samplenames <- colnames(readcounts)[4:ncol(readcounts)]
tumours <- intersect(samplenames, metadata[, tumour])
hosts <- intersect(samplenames, metadata[, host])
unknown <- setdiff(samplenames, union(tumours, hosts))

if (length(unknown) > 0) {
    logerror("Input readcounts contained unrecognised samples: %s",
             paste(unknown, collapse = ", "))
    stop()
} else {
    loginfo("No unrecognised samples found")
}

if (args$gccorrection != "") {
    loginfo("Doing GC correction")
    gc_content <- fread(args$gccorrection)
    if (ncol(gc_content) != 4) {
        stop(
            sprintf(
                "Expected to find 4 columns in GC interval table, but found %d",
                ncol(gc_content)))
    }
    colnames(gc_content) <- c("CHROM", "START", "END", "GC_CONTENT")
    logwarn(paste("Ensure that GC table START value is 1-based (i.e.",
            "intervals run from 10001-11000, not 10000-11000, for example)"))
    setkey(gc_content, CHROM, START, END)
    readcounts[gc_content, GC_CONTENT := i.GC_CONTENT]
    samplenames <- intersect(colnames(readcounts),
                             metadata[, c(tumour, host)])
    for (samplename in samplenames) {
        readcounts[, (samplename) := gatk_gc_correct(get(samplename),
                                                     GC_CONTENT)]
    }
}

loginfo("Processing read count matrices")
hosts_missing_tumours <-
    intersect(hosts, metadata[is.na(tumour) & !is.na(host), host])
male_hosts <-
    intersect(hosts, metadata[!is.na(host) & hostSex == "M", unique(host)])
female_hosts <-
    intersect(hosts, metadata[!is.na(host) & hostSex == "F", unique(host)])
sex_unknown <- setdiff(hosts, union(male_hosts, female_hosts))

if (length(hosts_missing_tumours) > 0) {
    logwarn("Found hosts with no matched tumour: %s",
            paste(hosts_missing_tumours, collapse = ", "))
}
if (length(sex_unknown) > 0) {
    logwarn("Found hosts with unknown sex: %s",
            paste(sex_unknown, collapse = ", "))
}

# Remove any noisy samples and/or those with polymorphic CNVs, or
# tumour contamination (based on metadata file)
hosts_to_exclude_from_panel <- vector("character")
if ("excludedFromPanel" %in% colnames(metadata)) {
    hosts_to_exclude_from_panel <- c(hosts_to_exclude_from_panel,
                                     metadata[excludedFromPanel == TRUE, host])
}
if ("tumourContaminated" %in% colnames(metadata)) {
    hosts_to_exclude_from_panel <- c(hosts_to_exclude_from_panel,
                                     metadata[tumourContaminated == TRUE, host])
}
hosts_to_exclude_from_panel <- unique(hosts_to_exclude_from_panel)
hosts <- setdiff(hosts, hosts_to_exclude_from_panel)

Hread_counts <- readcounts[, as.matrix(.SD), .SDcols = hosts]
Tread_counts <- readcounts[, as.matrix(.SD), .SDcols = tumours]
read_count_index <- readcounts[, .(CHROM, START, END)]

loginfo("Filtering out unsuitable bins")
loginfo("%s bins available", nrow(Hread_counts))

# Filters are TRUE where a data point should be retained,
# and FALSE where it should be removed (filtered out)
#
# All hosts have 10 or more reads in a bin
row_min_filter <- rowMins(Hread_counts) >= 10
# No host has more than 2000 reads in a bin
row_max_filter <- rowMaxs(Hread_counts) <= 2000
# Median read count must be above zero
row_medians_filter <- rowMedians(Hread_counts) >= 1
incomplete_bins_filter <- read_count_index[, (END - START + 1) %% 1000 == 0]

myfilter <- row_medians_filter &
    incomplete_bins_filter &
    row_min_filter &
    row_max_filter
loginfo("Filters removed %d bins, %d remaining", sum(!myfilter), sum(myfilter))

Hread_counts <- Hread_counts[myfilter, ]
Tread_counts <- Tread_counts[myfilter, ]
read_count_index <- read_count_index[myfilter, ]

stopifnot(nrow(Hread_counts) == sum(myfilter))
stopifnot(nrow(Tread_counts) == sum(myfilter))
stopifnot(nrow(read_count_index) == sum(myfilter))

loginfo("Saving filtered read counts to %s", args$outputpath)
fwrite(Hread_counts, file.path(args$outputpath,
                               "Hread_counts.csv.gz"))
fwrite(Tread_counts, file.path(args$outputpath,
                               "Tread_counts.csv.gz"))
fwrite(read_count_index, file.path(args$outputpath,
                               "read_count_index.csv.gz"))

loginfo("Creating the logR correctors")
all_hosts_corrector <- logRCorrector$new(Hread_counts)
female_hosts_corrector <-
    logRCorrector$new(Hread_counts[, colnames(Hread_counts) %in% female_hosts])
male_hosts_corrector <-
    logRCorrector$new(Hread_counts[, colnames(Hread_counts) %in% male_hosts])
loginfo("Saving logR correctors to %s", args$outputpath)
saveRDS(all_hosts_corrector, file.path(args$outputpath,
                                       "all_hosts_corrector.RDS"))
saveRDS(female_hosts_corrector, file.path(args$outputpath,
                                          "female_hosts_corrector.RDS"))
saveRDS(male_hosts_corrector, file.path(args$outputpath,
                                        "male_hosts_corrector.RDS"))

loginfo("Correcting the tumour logRs")
tumour_logr_standardized_all <-
    all_hosts_corrector$standardize(Tread_counts)
tumour_logr_standardized_male <-
    male_hosts_corrector$standardize(Tread_counts)
tumour_logr_standardized_female <-
    female_hosts_corrector$standardize(Tread_counts)

tumour_logr_denoised_all <- all_hosts_corrector$denoise(Tread_counts,
                                                        args$nsvs)
tumour_logr_denoised_male <- male_hosts_corrector$denoise(Tread_counts,
                                                          args$nsvs)
tumour_logr_denoised_female <- female_hosts_corrector$denoise(Tread_counts,
                                                              args$nsvs)

stopifnot(nrow(tumour_logr_denoised_all) == sum(myfilter))
stopifnot(nrow(tumour_logr_denoised_male) == sum(myfilter))
stopifnot(nrow(tumour_logr_denoised_female) == sum(myfilter))
stopifnot(all(tumours %in% colnames(tumour_logr_denoised_all)))
stopifnot(all(tumours %in% colnames(tumour_logr_denoised_male)))
stopifnot(all(tumours %in% colnames(tumour_logr_denoised_female)))

loginfo("Reassembling indexed count matrices")
Hread_counts <- cbind(read_count_index, Hread_counts)
Tread_counts <- cbind(read_count_index, Tread_counts)
tumour_logr_denoised_all <- cbind(read_count_index, tumour_logr_denoised_all)
tumour_logr_denoised_male <- cbind(read_count_index, tumour_logr_denoised_male)
tumour_logr_denoised_female <- cbind(read_count_index,
                                     tumour_logr_denoised_female)
setkey(Hread_counts, CHROM, START, END)
setkey(Tread_counts, CHROM, START, END)
setkey(tumour_logr_denoised_all, CHROM, START, END)
setkey(tumour_logr_denoised_male, CHROM, START, END)
setkey(tumour_logr_denoised_female, CHROM, START, END)

loginfo("Finding the host copy number")
# The inital part of this process builds a "row normaliser" to scale the
# host copy number to the correct integer.
# 1) Normalise the read counts by the column median
#    (median calculated over autosomes only)
read_count_index[, Index := .I]
autosome_index <- read_count_index[CHROM != "X" & CHROM != "Y", Index]
autosome_names <- read_count_index[CHROM != "X" & CHROM != "Y", unique(CHROM)]
allosome_index <- read_count_index[CHROM == "X" | CHROM == "Y", Index]
stopifnot(length(intersect(autosome_index, allosome_index)) == 0)
column_medians <- colMedians(Hread_counts[CHROM != "X" & CHROM != "Y",
                             as.matrix(.SD), .SDcols = hosts])
normalised_read_counts <-
    sweep(Hread_counts[, as.matrix(.SD), .SDcols = hosts], 2,
          column_medians, "/")

# 2) Use the row medians of the normalised columns to make the row normaliser,
#    for male and female samples
normalised_row_means <- rowMeans(normalised_read_counts)
normalised_row_means_male <-
    rowMeans(normalised_read_counts[,
             intersect(colnames(normalised_read_counts), male_hosts)])
normalised_row_means_female <-
    rowMeans(normalised_read_counts[,
             intersect(colnames(normalised_read_counts), female_hosts)])

male_row_normaliser <- normalised_row_means
male_row_normaliser[allosome_index] <-
    normalised_row_means_male[allosome_index]
male_row_normaliser[male_row_normaliser == 0] <- 1

female_row_normaliser <- normalised_row_means
female_row_normaliser[allosome_index] <-
    normalised_row_means_female[allosome_index]
female_row_normaliser[female_row_normaliser == 0] <- 1

# To find the copy number:
# 1) normalise the read counts: divide by the median count over the autosomes
# 2) divide by the sex-appropriate row normaliser
# 3) Take the log2 (fix infinities to an extreme finite number, e.g. 32)
# 4) Calculate the copy number using a purity & ploidy estimate of 1.0 and 2.0
# 5) Deal with sex chromosomes appropriately for the sample


loginfo("Writing outputs")
arrow_outdir <- file.path(args$outputpath, "dataset")
dir.create(arrow_outdir)

for (samplename_ in tumours) {
    print(samplename_)
    dt <- read_count_index[, .(CHROM, START, END)]
    matched_host <- metadata[tumour == samplename_, host]
    host_is_contaminated <- matched_host %in% hosts_to_exclude_from_panel
    host_sex <- metadata[tumour == samplename_, hostSex]
    dt[, samplename := samplename_]
    dt[, hostname := matched_host]
    dt[, T_bin_readcount := Tread_counts[, get(samplename_)]]

    if (is.na(matched_host)) {
        dt[, host_cn := NA]
        dt[, H_bin_readcount := NA_integer_]
    } else {
        dt[, H_bin_readcount := readcounts[myfilter, get(matched_host)]]

        # Calculate host copy number
        median_coverage <- dt[autosome_index, median(H_bin_readcount)]
        dt[, nm := H_bin_readcount / median_coverage]
        if (host_sex == "M") {
            dt[, nm := nm / male_row_normaliser]
        } else {
            dt[, nm := nm / female_row_normaliser]
        }
        dt[, logr := log2(nm)]
        dt[is.infinite(logr), logr := 32 * sign(logr)]

        if (host_sex == "M") {
            dt[, host_cn := cnpipe::calc_total_copynumber(logr, 1.0, 2.0,
                        ifelse(CHROM %in% autosome_names, 2, 1), 2.0)]
        } else {
            dt[, host_cn := cnpipe::calc_total_copynumber(logr, 1.0, 2.0,
                        ifelse(CHROM %in% c(autosome_names, "X"), 2, 0), 2.0)]
        }
        dt[, c("nm", "logr") := NULL]
    }

    # Annotate the sex-specific logR
    dt[, T_logr := tumour_logr_denoised_all[, get(samplename_)]]
    if (host_sex == "M") {
        dt[CHROM %in% c("X", "Y"),
           T_logr := tumour_logr_denoised_male[CHROM %in% c("X", "Y"),
                                               get(samplename_)]]
    } else {
        dt[CHROM %in% c("X", "Y"),
           T_logr := tumour_logr_denoised_female[CHROM %in% c("X", "Y"),
                                                 get(samplename_)]]
    }

    arrow::write_dataset(dt, arrow_outdir,
                         partitioning = "samplename",
                         format = "parquet")
    loginfo("Writing annotated host copy number for sample %s to %s",
            samplename_, args$outputpath)
}

loginfo("Done.")
