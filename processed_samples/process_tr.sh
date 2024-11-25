#!/bin/bash

for tuple in "bc2001 TRA" "bc2002 TRB" "bc2003 TRD" "bc2004 TRG"
#for tuple in "bc2001 TRA"
do
set -- $tuple
sample=$1
locus=$2
mkdir -p $(pwd)/results/$sample
mkdir -p $(pwd)/results/$sample/$locus
cd .
nextflow ../preprocess/main.nf --reads /mnt/f/clareo/easton_tr/m84248_241011_111920_s3.hifi_reads.$sample.fastq --outdir $(pwd)/results/$sample/$locus --sample_name $sample --locus $locus --constant_region TR
cd .
nextflow ../annotate_tr/main.nf --reads $(pwd)/results/$sample/$locus/reads/${sample}_atleast-2.fasta --outdir $(pwd)/results/$sample/$locus --sample_name $sample --locus $locus
done

cd ..
docker run -i --cpus 60.0 -v $(pwd):/data -w /data williamlees/flairr_dsl2:latest /data/processed_samples/run_summaries.sh