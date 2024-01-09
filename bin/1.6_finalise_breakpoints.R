library(logging)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("inputpath", help = "Path to processed logR data")
parser$add_argument("metadata", help = "Path to metadata file [Delimited format]")
# parser$add_argument("clone", help = "Clone for which to combine segment IDs")
parser$add_argument("savepath", help="Path to write final output. Will be created if not found")
args <- parser$parse_args()

# Example
# args <- list(
    # inputpath = "../../CTVT/output/2022-06-15_CopyWright_1.3_filter_breakpoints",
    # metadata="../../CTVT/data/dog_metadata.csv",
    # clone="CTVT",
    # savepath = "../../CTVT/output/2022-06-15_CopyWright_1.4_finalise_breakpoints"
# )


for (argname in names(args)) {
    logdebug("%s = %s", argname, args[[argname]])
}

if (!dir.exists(args$inputpath)) {
    stop(sprintf("Could not find input path %s", args$inputpath))
}
if (!dir.exists(file.path(args$inputpath, "dataset"))) {
    stop(sprintf("Could not find dataset path %s",
                 file.path(args$inputpath, "dataset")))
}
if (!dir.exists(args$savepath)) {
    dir.create(args$savepath, recursive=TRUE)
}
if (!dir.exists(args$savepath)) {
    stop(sprintf("Could not create save path %s", args$savepath))
}
if (!file.exists(args$metadata)) {
    stop(sprintf("Could not find metadata file %s", args$metadata))
}

suppressPackageStartupMessages({
    library(data.table)
    library(arrow)
    library(dplyr)
    library(progress)
})

loginfo("Loading data")
metadata <- fread(args$metadata, na.strings="")
# if (!(args$clone %in% metadata$clone)) {
#     stop(sprintf("No data available for clone %s", args$clone))
# }


segmentation_files <- list.files(args$inputpath, "*.csv")

all_breakpoints <- rbindlist(lapply(file.path(args$inputpath,
                                              segmentation_files), fread))
unified_breakpoints <- unique(all_breakpoints,
                              by = c("CHROM", "START", "clone"))[order(clone,
                                                                       CHROM,
                                                                       START)]

data_dir <- file.path(args$inputpath, "dataset")
ds <- arrow::open_dataset(data_dir)

samplenames <- metadata[, samplename]

pb <- progress_bar$new(format="Finalising calls [:bar] :percent :eta",
                       total = length(samplenames),
                       width = 75,
                       clear = FALSE)

loginfo("Finalising segmentation calls")
for (samplename_ in samplenames) {
    ds %>% filter(samplename == samplename_) %>% collect() -> dt
    setDT(dt)
    setcolorder(dt, c("CHROM", "START", "END", "samplename"))
    setorder(dt, CHROM, START, END, samplename)
    clone_ <- metadata[samplename == samplename_, clone]
    dt[, dataset_segment_id := 0]
    dt[, is_bk := 0]
    dt[unified_breakpoints[clone==clone_], is_bk := 1,
       on = c("CHROM", "START")]
    dt[, dataset_segment_id := cumsum(is_bk), by = CHROM]
    dt[dataset_segment_id == 0, dataset_segment_id := 1]
    dt[, pipeline_cn := round(median(ifelse(excluded_from_segmentation, NA,
                                            total_cn),
                                     na.rm = TRUE)),
       by = .(CHROM, dataset_segment_id)]

    outdir <- file.path(args$savepath, "dataset")
    setorder(dt, CHROM, START, END, samplename)
    arrow::write_dataset(dt[, .(CHROM, START, END, samplename, hostname,
                                T_bin_readcount, H_bin_readcount,
                                excluded_from_segmentation,
                                host_cn, expected_host_cn, T_logr, total_cn,
                                sample_segment_id=segment_id,
                                sample_segment_cn=segment_median,
                                dataset_segment_id, pipeline_cn)],
                         outdir, partitioning = "samplename")
    pb$tick()
}
