#!/bin/bash

# nextflow ../preprocess/main.nf --reads /mnt/f/clareo/easton_tr/test_IGK.fastq --outdir $(pwd)/results/IGK --sample_name test_IGK --locus IGK --constant_region IGK
cd .
nextflow ../annotate_tr/main.nf --reads /mnt/f/clareo/easton_tr/m84248_241011_111920_s3.hifi_reads.bc2001.fasta --outdir $(pwd)/results/TRA --sample_name bc2001 --locus TRA
cd .
nextflow ../annotate_tr/main.nf --reads /mnt/f/clareo/easton_tr/m84248_241011_111920_s3.hifi_reads.bc2002.fasta --outdir $(pwd)/results/TRB --sample_name bc2002 --locus TRB
cd .
nextflow ../annotate_tr/main.nf --reads /mnt/f/clareo/easton_tr/m84248_241011_111920_s3.hifi_reads.bc2003.fasta --outdir $(pwd)/results/TRD --sample_name bc2003 --locus TRD
cd .
nextflow ../annotate_tr/main.nf --reads /mnt/f/clareo/easton_tr/m84248_241011_111920_s3.hifi_reads.bc2004.fasta --outdir $(pwd)/results/TRG --sample_name bc2004 --locus TRG
