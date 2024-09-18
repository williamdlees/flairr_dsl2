#!/bin/bash

# Example script to process multiple repertoires
# The repertoires, and personalised reference sets, where specified, must be given as full (not realtive) pathnames
# The script creates a directory underneath processed_samples for each repertoire, and runs the pipeline in that directory

# Example using 'standard' reference sets. The locations are specified in processed_samples/FLAIRRSeq.config

for sample in 1001
do
mkdir $sample-kinnex
cd $sample-kinnex
nextflow ../../preprocess/main.nf --reads /mnt/f/clareo/kinnex_samples/m84248_240915_090555_s3.hifi_reads.RACE_TSO_5p--bc*$sample*.hifi_reads.fastq --outdir $(pwd)/results
# cd is needed to avoid a Java bug - probably wsl specific
cd ..
cd $sample-kinnex
nextflow ../../annotate/main.nf --reads $(pwd)/results/reads/*collapse-unique_atleast-2.fasta --outdir $(pwd)/results
cd ..
docker run -v "$(pwd)/..":/scratch williamlees/flairr_dsl2 /scratch/processed_samples/presto_r.sh $sample-kinnex
done

cd ..
docker run -i --cpus 60.0 -v $(pwd):/data -w /data my_immcantation_suite:4.4.0 /data/processed_samples/run_summaries.sh
