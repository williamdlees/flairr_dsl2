#!/bin/bash
cd /scratch/processed_samples
R -e "rmarkdown::render('../presto_r/FLAIRR.Rmd', params=list(data='../processed_samples/"$1"/results/reports', sample='"$1"'), output_file='../processed_samples/"$1"/"$1"_FLAIRR_presto.pdf')"
