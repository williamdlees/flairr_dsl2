#!/bin/bash

pip install toml
pip install plotly

cd processed_samples
python3 ../python/pipeline_counts.py flairr_logs.toml pipeline_counts.csv
python3 ../python/barcode_alignments.py flairr_logs.toml barcode_alignments.html
