
nextflow ../preprocess/main.nf --reads /mnt/f/clareo/easton_samples/IGL/demultiplex.RACE_TSO_5p--bc1001.bam.fasta.fastq --outdir $(pwd)/results/IGL --sample_name test_IGL --locus IGL --constant_region IGL
cd .
nextflow ../annotate/main.nf --reads /mnt/f/clareo/flairr_dsl2/processed_samples/results/IGL/reads/test_IGL_atleast-2.fasta --outdir $(pwd)/results/IGL --sample_name IGL --locus IGL
cd .
#nextflow ../annotate_with_inference/main.nf --sample_name test --reads /mnt/f/clareo/flairr_dsl2/processed_samples/results/reads/test_atleast-2.fasta --outdir $(pwd)/results_with_inference
cd .
