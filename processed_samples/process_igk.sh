
nextflow ../preprocess/main.nf --reads /mnt/f/clareo/easton_samples/IGK/m64407e_230922_095324.hifi_reads.fastq --outdir $(pwd)/results --sample_name test_IGK --locus IGK --constant_region IGK
cd .
nextflow ../annotate/main.nf --reads /mnt/f/clareo/flairr_dsl2/processed_samples/results/reads/test_IGK_atleast-2.fasta --outdir $(pwd)/results --sample_name test_IGK --locus IGK --haplotype_genes IGKJ2
cd .
#nextflow ../annotate_with_inference/main.nf --sample_name test --reads /mnt/f/clareo/flairr_dsl2/processed_samples/results/reads/test_atleast-2.fasta --outdir $(pwd)/results_with_inference
cd .
