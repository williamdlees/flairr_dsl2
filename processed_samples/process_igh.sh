#!/bin/bash

# Example script to process multiple repertoires
# The repertoires, and personalised reference sets, where specified, must be given as full (not realtive) pathnames
# The script creates a results directory underneath processed_samples for each repertoire, and runs the pipeline in that directory

# Example using 'standard' reference sets. The locations are specified in processed_samples/FLAIRRSeq.config

#for sample in 1001 1002 1003 1004 1005 1006 1008 1012
for sample in 1001
do
#nextflow ../preprocess/main.nf --reads /mnt/f/clareo/easton_samples/*667*$sample*.fastq --outdir $(pwd)/results/IGH --sample_name $sample-667
# cd is needed to avoid a Java bug - probably wsl specific
cd .
nextflow ../annotate/main.nf --reads $(pwd)/results/IGH/reads/${sample}-667_atleast-2.fasta --outdir $(pwd)/results/IGH  --sample_name $sample-667
cd .
done

# Example using 'personalised' V reference sets. The specification on the command line overrides processed_samples/FLAIRRSeq.config. D and J could also be specified.
# The v_germline file must be IMGT gapped

for sample in 1001 1003
do
#nextflow ../preprocess/main.nf --reads /mnt/f/clareo/easton_samples/*986*$sample*.fastq --outdir $(pwd)/results/IGH --sample_name $sample-986
# cd is needed to avoid a Java bug - probably wsl specific
cd .
#nextflow ../annotate/main.nf --reads $(pwd)/results/IGH/reads/${sample}-986_atleast-2.fasta --outdir $(pwd)/results/IGH --sample_name $sample-986
done

cd ..
#docker run -i --cpus 60.0 -v $(pwd):/data -w /data my_immcantation_suite:4.4.0 /data/processed_samples/run_summaries.sh
