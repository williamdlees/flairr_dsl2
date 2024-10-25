#!/bin/bash

cd /data/processed_samples
python3 ../python/pipeline_counts.py flairr_logs.toml results/pipeline_counts.csv
python3 ../python/barcode_alignments.py flairr_logs.toml results/barcode_alignments.html
