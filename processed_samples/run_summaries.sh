#!/bin/bash

# Start in the mounted data directory (processed_samples)
cd /data

# Run each summary-generating script with absolute paths
echo "🔄 Running pipeline_counts.py..."
python3 /nf_root/python/pipeline_counts.py /nf_root/preprocess/flairr_logs.toml /data/results/pipeline_counts.csv

echo "🔄 Running barcode_alignments.py..."
python3 /nf_root/python/barcode_alignments.py /nf_root/preprocess/flairr_logs.toml /data/results/barcode_alignments.html

echo "🔄 Running project_allele_table.py..."
python3 /nf_root/python/project_allele_table.py flairr_test /data/results /data/results/flairr_allele_usage_table.csv

echo "🔄 Running make_project_alignments_file.py..."
cd /data/results
#python3 /nf_root/python/make_project_alignments_file.py

echo "✅ All summaries complete."
