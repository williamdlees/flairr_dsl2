## Demultiplexing PacBio FLAIRR-seq data

This note explains how to demultiplex the single BAM file produced by the PacBio sequencer into separate FASTQ files for each sample. The demultiplexing is done in two steps: first, the reads are demultiplexed into pools based on the pool barcodes, and then each pool is demultiplexed into individual samples based on the sample barcodes.

Preparation:
- make a directory under /mnt/efs, for example `./pools`. In that directory:
- In that directory, create a file for each pool specifying the expected barcodes in each read and the corresponding sample name. Follow the format of the example at
[pool_A_biosample.csv](pool_A_biosample.csv). 
- Create a subdirectory for the sequencing run, for example `sequencing`. 
- Create a subdirectory for each pool. You can use any names you like but I tend to use the pool barcode names, for example bc2048, bc2049 etc. 

The first step is to demultiplex into pools. This is done with PacBio lima. cd to the sequencing run directory and run the command:

```bash
lima m21114_260402_212026.hifi_reads.bam smrt_adapters.fasta output.demux.bam --same --split-bam-named --min-score 95
```

`m21114_260402_212026.hifi_reads.bam` is the path to the file provided by the sequencer. Note that `min_score` is set to 95, which is higher than the default. 
For simplicity, you may wish to rename the output files to a simple format that includes the pool barcode, for example `bc2041.bam`, `bc2042.bam` etc. This is not essential but it makes it easier to keep track of the files. If you do that, remember to rename the `.pbi` files also.
[smrt_adapters.fasta](smrt_adapters.fasta) is the file containing the PacBio adapter sequences.

When lima has completed, cd to the first of the pool directories you created and run this command:

```bash
python python/simple_demux.py \
 --input ../sequencing/bc2041.bam \
 --barcodes barcodes.fasta \
 --biosample ../bc2041.csv  \
 --chunk-size 100000
```

`--input` specifies the bam file for the pool, created in the last step. `--biosample` specifies the file you created that maps the barcodes to sample names for this pool. `--barcodes` specifies the file containing the expected barcode sequences. The [example file](idt_barcodes.fasta) linked here contains IDT barcodes, currently used in the FLAIRR protocol.

The demultiplexing tool will remove from each read the barcode and TSO sequences listed in the barcode file. The 3-prime barcode sequences should therefore not include the primer sequence, as this is required to be present in the read for later processing. 

The command will create demultiplexed FASTQ files in the directory. At this point you should check that the output from the tool (`written to demux_report.txt`) is reasonable and the sizes of the files are reasonable. You can also review `audit_log.csv`
which gives a read-by-read summary of what the demuxer found.

Repeat this step for each pool.