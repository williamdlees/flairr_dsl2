# FLAIRR DSL2 Pipeline

A Nextflow-based pipeline for preprocessing and annotating immunoglobulin and T-cell FLAIRR-seq receptor repertoires, designed for high-performance computing environments using Slurm workload management.

## Overview

The FLAIRR DSL2 pipeline provides a comprehensive workflow for analyzing adaptive immune receptor repertoires from high-throughput sequencing data produced by the FLAIRR-seq protocol. The pipeline supports both immunoglobulin (IGH, IGK, IGL) and T-cell receptor (TRA, TRB, TRD, TRG) loci, offering:

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

The main entry point is the `process_slurm.sh` script, which orchestrates batch processing of repertoire samples. 

```bash
# Preprocess samples
./process_slurm.sh preprocess input_samples.txt IGH 5 docker

# Annotate preprocessed samples  
./process_slurm.sh annotate sample_list.txt TRB 10 singularity -p bigmem
```

**Documentation**

Complete usage information for process_slurm.sh is available in the man page:
```bash
man -l processed_samples/process_slurm.1
```

You can also view the processed man page as plain text: [process_slurm_man.txt](processed_samples/process_slurm_man.txt)

The output (results from processing and logs of the processing run) will be stored in the current working directory. The script itself is located in the `processed_samples/` directory of the repository. You can either call it from there or cd to some other preferred working directory and call it with a full or relative path. Best practice, if you intend to run multiple projects, is to make a new empty directory for each project and call the script from there.

**Key Features:**
- **Batch Processing**: Submits individual Slurm jobs for each sample
- **Job Throttling**: Manages concurrent job limits to respect cluster policies
- **Flexible Configuration**: Supports custom Slurm partitions and CPU allocation
- **Container Runtime**: Works with both Docker and Singularity
- **Parameter Passing**: Forwards additional Nextflow parameters seamlessly

**Input Format:**
Tab-separated file with sample names and FASTQ paths:
```
sample_01	/path/to/sample_01.fastq
sample_02	/path/to/sample_02.fastq
```

**Output Directories:**
- `./results/`: Pipeline results organized by sample and locus
- `./slog/`: Slurm job logs for monitoring and debugging

### 2. Project Summaries with `singularity_summaries.sh`

After all samples have been processed and annotated, generate project-wide summaries:

```bash
./singularity_summaries.sh
```

For Docker usage, the script `run_summaries.sh` can be run from within the flairr_dsl2 container, with suitable mount points to the results directory and the root of the repository.

This script should be run from the same directory that you ran process_slurm.sh. (adjust the path to singularity_summaries.sh appropriately).

This script produces:
- **Pipeline Statistics**: Processing counts and success rates
- **Barcode Alignments**: Quality control visualizations  
- **Allele Usage Tables**: Gene usage patterns across the project
- **Project Alignments**: Consolidated alignment files

## Container Setup

### Docker
The pipeline uses pre-built Docker containers. Ensure Docker is available and properly configured on your system. Containers will be pulled automatically when the pipeline is executed.

### Singularity
For Singularity usage, convert the referenced Docker containers to `.sif` files (see [shared_configs/process.config](shared_configs/process.config) for a list of containers). 

**Note:** the scripts expect these files to be stored in the directory `~/sifs/`: if you need to keep them somewhere else, please link `~/sifs/` to that location. 

```bash
# Example conversion
singularity build flairr_dsl2_latest.sif docker:~/sifs/flairr_dsl2:latest
```

Refer to the [Singularity documentation](https://docs.sylabs.io/guides/latest/user-guide/) for detailed conversion procedures and best practices.


## Quick Start

1. **Clone the Repository**: Use Git to clone the repository to your local machine.
2. **Build the .sif files**: If using Singularity, convert the Docker images to Singularity .sif files.
3. **Make an empty directory**: cd to it and call the shell scripts from this location.
4. **Prepare Input**: Create a tab-separated file listing your samples and FASTQ paths.
5. **Reference directory**: Make a reference directory containing the necessary files (See below).
5. **Preprocess**: Run `process_slurm.sh preprocess` to preprocess your data.
6. **Annotate**: Run `process_slurm.sh annotate` to perform V(D)J assignment.
7. **Summarize**: Execute `singularity_summaries.sh` to generate project reports.

## Requirements

- **Slurm Workload Manager**: For job scheduling and resource management
- **Nextflow**: Workflow execution engine
- **Container Runtime**: Docker or Singularity
- **Python 3**: For summary script execution

## Reference Directory

For each locus, the reference directory should contain the following files:
- Homo_sapiens_\<locus\>V.fasta  
- Homo_sapiens_\<locus\>D.fasta (if applicable)
- Homo_sapiens_\<locus\>J.fasta
- Homo_sapiens_\<locus\>VDJ.fasta (concatentation of the above three files)
- Homo_sapiens_\<locus\>V_gapped.fasta - IMGT-aligned V germline sequences
- Homo_sapiens_\<locus\>C.fasta
- Homo_sapiens_\<locus\>.aux - IGBLAST-format aux file
- Homo_sapiens_\<locus\>.ndm - IGBLAST-format ndm file

## Support

For detailed usage instructions and troubleshooting:
- View the manual: `man -l processed_samples/process_slurm.1`
- Check Slurm logs in `./slog/` directory
- Review Nextflow logs within individual sample result directories

## FAQs

## Running slurm on a local system

You can install and run slurm on a single machine. This is a useful and practical way to run jobs locally, with full management and tracking across multiple samples. Slurm is supported on most linux distributions and also runs on Windows via WSL2. Check the slurm documentation for MacOS support - there were issues at the time of writing this README. Notes on setting up slurm on a single machine can be found [here](https://drtailor.medium.com/how-to-setup-slurm-on-ubuntu-20-04-for-single-node-work-scheduling-6cc909574365).

## Running the pipeline under Windows WSL2

The installation of `muscle` in the Immcantation container uses static libraries not available under WSL2. Unless these are fixed, preprocessing will freeze at the AlignSets step. To resolve this, log in to the container and rebuild `muscle`, without the static flag.

```bash
docker run -it --cpus 1 immcantation/suite:4.4.0 /bin/bash
mkdir muscledir
cd muscledir
wget https://www.drive5.com/muscle/muscle_src_3.8.1551.tar.gz
tar xzvf muscle_src_3.8.1551.tar.gz
make
# the make will fail. Now run the last step again (the linker command which failed) without the -static option

# test 
muscle --help

# move the new muscle binary to /usr/local/bin to override the existing one
mv muscle /usr/local/bin
cd ..
rm -rf muscledir
```

At this point, use another command window to find the container id and commit the changes to a new image.

```bash
docker ps
docker commit <container_id> immcantation/suite:4.4.0_wsl2
```

Now modify `./shared_configs/process.config` to use `immcantation/suite:4.4.0_wsl2` as the container image for the immcantation processes. You should now be able to run the pipeline under WSL2.

## Citation

When using this pipeline, please cite the following publication.

**FLAIRR-seq Method:**
Ford EE, Tieri D, Rodriguez OL, Francoeur NJ, Soto J, Kos JT, Peres A, Gibson WS, Silver CA, Deikus G, et al. FLAIRR-Seq: A Method for Single-Molecule Resolution of Near Full-Length Antibody H Chain Repertoires. The Journal of Immunology (2023) 210:1607–1619. doi: [10.4049/jimmunol.2200825](https://doi.org/10.4049/jimmunol.2200825)

