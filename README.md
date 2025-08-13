# FLAIRR DSL2 Pipeline

A Nextflow-based pipeline for preprocessing and annotating immunoglobulin and T-cell receptor repertoires, designed for high-performance computing environments using Slurm workload management.

## Overview

The FLAIRR DSL2 pipeline provides a comprehensive workflow for analyzing adaptive immune receptor repertoires from high-throughput sequencing data. The pipeline supports both immunoglobulin (IGH, IGK, IGL) and T-cell receptor (TRA, TRB, TRD, TRG) loci, offering:

- **Preprocessing**: Quality control, barcode processing, and read preparation
- **Annotation**: V(D)J gene assignment, CDR3 identification, and repertoire analysis
- **Summary Generation**: Project-wide statistics and visualizations

## Architecture

The pipeline is designed specifically for **Slurm-managed HPC environments** and supports containerized execution using either:
- **Docker** containers
- **Singularity** containers (.sif files)

All computational processes run within containers to ensure reproducibility and environment consistency across different systems.

## Core Workflow

### 1. Sample Processing with `process_slurm.sh`

The main entry point is the `process_slurm.sh` script, which orchestrates batch processing of repertoire samples:

```bash
# Preprocess samples
./process_slurm.sh preprocess input_samples.txt IGH 5 docker

# Annotate preprocessed samples  
./process_slurm.sh annotate sample_list.txt TRB 10 singularity -p bigmem
```

**Key Features:**
- **Batch Processing**: Submits individual Slurm jobs for each sample
- **Job Throttling**: Manages concurrent job limits to respect cluster policies
- **Flexible Configuration**: Supports custom Slurm partitions and CPU allocation
- **Container Runtime**: Works with both Docker and Singularity
- **Parameter Passing**: Forwards additional Nextflow parameters seamlessly

**Input Format:**
Tab-separated file with sample names and FASTQ paths:
```
sample_01	/path/to/sample_01.fastq.gz
sample_02	/path/to/sample_02.fastq.gz
```

**Output Directories:**
- `./results/`: Pipeline results organized by sample and locus
- `./slog/`: Slurm job logs for monitoring and debugging

**Documentation:**
Complete usage information is available in the man page:
```bash
man -l processed_samples/process_slurm.1
```

You can also view the processed man page as plain text: [process_slurm_man.txt](processed_samples/process_slurm_man.txt)

### 2. Project Summaries with `singularity_summaries.sh`

After all samples have been processed and annotated, generate project-wide summaries:

```bash
./processed_samples/singularity_summaries.sh
```

This script produces:
- **Pipeline Statistics**: Processing counts and success rates
- **Barcode Alignments**: Quality control visualizations  
- **Allele Usage Tables**: Gene usage patterns across the project
- **Project Alignments**: Consolidated alignment files

## Container Setup

### Docker
The pipeline uses pre-built Docker containers. Ensure Docker is available and properly configured on your system.

### Singularity
For Singularity usage, convert the referenced Docker containers to `.sif` files. Note that the scripts expect these files to be in `~/sifs/`: if you need t store them somewhere else please link ~/sifs to that location. 

```bash
# Example conversion
singularity build flairr_dsl2_latest.sif docker:~/sifs/flairr_dsl2:latest
```

Refer to the [Singularity documentation](https://docs.sylabs.io/guides/latest/user-guide/) for detailed conversion procedures and best practices.

**Note:** Update the singularity image path in `singularity_summaries.sh` to match your local `.sif` file location.

## Quick Start

1. **Prepare Input**: Create a tab-separated file listing your samples and FASTQ paths
2. **Preprocess**: Run `process_slurm.sh preprocess` to prepare your data
3. **Annotate**: Run `process_slurm.sh annotate` to perform V(D)J assignment
4. **Summarize**: Execute `singularity_summaries.sh` to generate project reports

## Requirements

- **Slurm Workload Manager**: For job scheduling and resource management
- **Nextflow**: Workflow execution engine
- **Container Runtime**: Docker or Singularity
- **Python 3**: For summary script execution

## Support

For detailed usage instructions and troubleshooting:
- View the manual: `man -l processed_samples/process_slurm.1`
- Check Slurm logs in `./slog/` directory
- Review Nextflow logs within individual sample result directories

## Citation

When using this pipeline, please cite the appropriate publications for the underlying tools and methodologies implemented in the FLAIRR framework.

