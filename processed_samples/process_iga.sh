#!/bin/bash

for sample in 5104
do
nextflow ../preprocess/main.nf --reads /mnt/f/clareo/easton_IgA/$sample*.fastq --outdir $(pwd)/results/IGH --sample_name $sample-IgA --constant_region AGM --align_sets.nproc 1 -resume
# cd is needed to avoid a Java bug - probably wsl specific
cd .
# nextflow ../annotate/main.nf --reads $(pwd)/results/IGH/reads/${sample}-IgA_atleast-2.fasta --outdir $(pwd)/results/IGH  --sample_name $sample-IgA
cd .
done

