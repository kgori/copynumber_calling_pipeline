copynumber_calling_pipeline
==

## Summary
Nextflow pipeline for calling copynumber variants from transmissible cancer samples.

Used in the nuclear horizontal transfer in CTVT project.

## Dependencies
I've not made a container for this one. It requires that the following R packages are installed:
 - argparse
 - arrow
 - [cnpipe](https://github.com/kgori/cnpipe.git)
 - data.table
 - dplyr
 - logging
 - magrittr
 - matrixStats
 - naturalsort
 - progress
 - randomForest
 - [segmentation](https://github.com/kgori/segmentation.git)

`cnpipe` and `segmentation` are my packages. They can be installed using `devtools`/`remotes::install_github`.
The package repos have installation details.
