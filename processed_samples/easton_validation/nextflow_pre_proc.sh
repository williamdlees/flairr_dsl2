#!/bin/bash

# Arguments
input_file=$1

# Loop through the TSV
while IFS=$'\t' read -r sampleId pathToReads; do
    mkdir -p ./results/${sampleId}/IGH
    # Check the queue count
    count=$(squeue -p bioinfo | grep $USER | wc -l)

    # Limit job submission to 10 jobs at a time
    while [ "$count" -gt 10 ]; do
        sleep 1s
        count=$(squeue -p bioinfo | grep $USER | wc -l)
    done

    # Submit the job
    sbatch -p bioinfo \
       --oversubscribe \
       -J "${sampleId}_preproc" \
       -o "${sampleId}_preproc.slog" \
       --cpus-per-task=12 \
       --wrap="nextflow run -offline /home/w0lees01/flairr_dsl2/preprocess/main.nf \
            --sample_name '${sampleId}' \
            --reads '${pathToReads}' \
            --outdir './results/${sampleId}/IGH'"  
done < "$input_file"
