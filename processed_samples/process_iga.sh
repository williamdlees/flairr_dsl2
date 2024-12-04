#!/bin/bash

for sample in 1001_m84248_241112_211045_s2.hifi_reads.RACE_TSO_5p--bc1004.bam_trunc
do
cd .
nextflow ../preprocess/main.nf --reads /mnt/f/clareo/easton_IgA/$sample*.fastq --outdir $(pwd)/results/$sample/IGH --sample_name $sample-IgA --constant_region AGM --nproc 40 --align_sets.nproc 8 -resume
# cd is needed to avoid a Java bug - probably wsl specific
cd .
nextflow ../annotate/main.nf --reads $(pwd)/results/$sample/IGH/reads/${sample}-IgA_atleast-2.fasta --outdir $(pwd)/results/$sample/IGH  --sample_name $sample-IgA
cd .
done

