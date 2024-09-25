
nextflow ../preprocess/main.nf --reads /mnt/f/clareo/easton_short_samples/*986*1003*.fastq --outdir $(pwd)/results --sample_name test
cd .
nextflow ../annotate_with_inference/main.nf --sample_name test --reads /mnt/f/clareo/flairr_dsl2/processed_samples/results/reads/test_atleast-2.fasta --outdir $(pwd)/results_with_inference
cd .
