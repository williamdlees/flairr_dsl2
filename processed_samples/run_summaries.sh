#!/bin/bash

cd /data/processed_samples
python3 ../python/pipeline_counts.py flairr_logs.toml results/pipeline_counts.csv
python3 ../python/barcode_alignments.py flairr_logs.toml results/barcode_alignments.html
python3 ../python/project_allele_table.py flairr_test, results, results/flairr_allele_usage_table.csv
