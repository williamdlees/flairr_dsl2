#!/bin/bash
# process_slurm.sh
# Usage: ./process_slurm.sh <command> <input_fofn> <locus> <max_jobs> <container_runtime> [additional_nextflow_parameters...]
# <command> must be preprocess or annotate
# <input_fofn> should consist of lines of the form <sample_name>\t<fastq_path_name>
# <container_runtime> must be docker or singularity
# [additional_nextflow_parameters...] will be passed directly to nextflow
# slurm log files are written to the subdirectory slog

# Function to display help message
show_help() {
  cat << EOF
USAGE:
  ./process_slurm.sh <command> <input_fofn> <locus> <max_jobs> <container_runtime> <germline_ref_dir> [options] [additional_nextflow_parameters...]

DESCRIPTION:
  This script processes samples in batch through Nextflow pipelines using Slurm.
  
REQUIRED ARGUMENTS:
  <command>             Must be either 'preprocess' or 'annotate'
  <input_fofn>          Tab-separated file with <sample_name> and <fastq_path_name> on each line
  <locus>               Target locus: IGH, IGK, IGL, TRA, TRB, TRD, or TRG
  <max_jobs>            Maximum number of jobs to run in parallel
  <container_runtime>   Container runtime to use: 'docker' or 'singularity'
  <germline_ref_dir>    Path to the germline reference directory

OPTIONS:
  -p <partition>        Slurm partition to use (default: bioinfo)
  --process.cpus <n>    Number of CPUs to request per job (default: 12)
                        This is also passed to Nextflow for process configuration
  --echo                Echo commands instead of executing them
  -help, --help         Display this help message

EXAMPLES:
  ./process_slurm.sh preprocess input_samples.txt IGH 5 docker /path/to/germline/refs
  ./process_slurm.sh annotate results.txt TRB 10 singularity /path/to/germline/refs -p bigmem --process.cpus 16
  ./process_slurm.sh preprocess input_samples.txt IGH 5 docker /path/to/germline/refs --echo

OUTPUT:
  - Results are stored in the './results/<sample>/' directory
  - Slurm logs are written to the './slog/' directory
EOF
  exit 0
}

# Check for help flags first
for arg in "$@"; do
  if [[ "$arg" == "-help" || "$arg" == "--help" ]]; then
    show_help
  fi
done

set -euo pipefail

# setup paths to nextflow scripts
SCRIPT_PATH="${BASH_SOURCE[0]}"
# Convert to absolute path
SCRIPT_PATH="$(readlink -f "$SCRIPT_PATH")"
SCRIPT_DIR="$(dirname "$SCRIPT_PATH")"
NF_ROOT="${SCRIPT_DIR%/*}"
readonly NXF_PP_SCRIPT="$NF_ROOT/preprocess/main.nf"
readonly NXF_IG_SCRIPT="$NF_ROOT/annotate/main.nf"
readonly NXF_TR_SCRIPT="$NF_ROOT/annotate_tr/main.nf"

mkdir -p slog
mkdir -p results

# Default partition
partition="bioinfo"
# Default CPU count
cpus_per_task=12
# Default echo mode
echo_only=false

# Process parameters
USAGE="Usage: $0 <command> <input_tsv> <locus> <max_jobs> <container_runtime (docker or singularity)> <germline_ref_dir> [additional_nextflow_parameters...]. Use $0 --help for help."
command="${1:?{$USAGE}"
input_file="${2:?{$USAGE}"
locus="${3:?$USAGE}"
max_jobs="${4:?$USAGE}"
runtime="${5:?$USAGE}"
germline_ref_dir="${6:?$USAGE}"

# Shift the first 6 mandatory parameters
shift 6

# Check for partition option and process.cpus
additional_params=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        -p)
            partition="$2"
            shift 2
            ;;
        --process.cpus)
            if [[ $2 =~ ^[0-9]+$ ]]; then
                cpus_per_task="$2"
                additional_params+="$1 $2 "
                shift 2
            else
                echo "Error: --process.cpus requires a numeric argument"
                exit 1
            fi
            ;;
        --echo)
            echo_only=true
            shift
            ;;
        *)
            additional_params+="$1 "
            shift
            ;;
    esac
done

# Check if command has a valid value
valid_commands=("preprocess" "annotate")
valid_command=false
for valid in "${valid_commands[@]}"; do
    if [ "$command" = "$valid" ]; then
        valid_command=true
        break
    fi
done


if [ "$valid_command" = false ]; then
    echo "Error: Invalid locus '$command'. Must be one of: preprocess, annotate."
    exit 1
fi

# Check if input file exists
if [ ! -f "$input_file" ]; then
    echo "Error: Input file '$input_file' does not exist."
    exit 1
fi

# Check if locus has a valid value
valid_loci=("IGH" "IGK" "IGL" "TRA" "TRB" "TRD" "TRG")
valid_locus=false
for valid in "${valid_loci[@]}"; do
    if [ "$locus" = "$valid" ]; then
        valid_locus=true
        break
    fi
done

if [ "$valid_locus" = false ]; then
    echo "Error: Invalid locus '$locus'. Must be one of: IGH, IGK, IGL, TRA, TRB, TRD, TRG."
    exit 1
fi

# Check if germline_ref_dir exists
if [ ! -d "$germline_ref_dir" ]; then
    echo "Error: Germline reference directory '$germline_ref_dir' does not exist."
    exit 1
fi

# Check if max_jobs is an integer between 1 and 25
if ! [[ "$max_jobs" =~ ^[0-9]+$ ]]; then
    echo "Error: max_jobs '$max_jobs' must be an integer."
    exit 1
fi

if [[ "$command" == "annotate" ]]; then
  if [[ "$locus" == *"TR"* ]]; then
      NXF_SCRIPT=$NXF_TR_SCRIPT
  else
      NXF_SCRIPT=$NXF_IG_SCRIPT
      pathToReads=
  fi
else
  NXF_SCRIPT=$NXF_PP_SCRIPT
fi

# Check if runtime has a valid value
valid_runtimes=("docker" "singularity")
valid_runtime=false
for valid in "${valid_runtimes[@]}"; do
    if [ "$runtime" = "$valid" ]; then
        valid_runtime=true
        break
    fi
done

if [ "$valid_runtime" = false ]; then
    echo "Error: Invalid container runtime '$runtime'. Must be one of: docker, singularity."
    exit 1
fi

# check the runtime is available

if command -v $runtime > /dev/null 2>&1; then
    echo " "
else
    echo "$runtime is not available. Maybe you need to load a module?"
    exit 1
fi


# Check that script file exists
if [ ! -f "$NXF_SCRIPT" ]; then
    echo "Error: Script file '$NXF_SCRIPT' does not exist."
    exit 1
fi


while IFS=$'\t' read -r sample pathToReads; do
  # skip empty or malformed lines
  [[ -z "$sample" ]] && continue

  # sanitize sample name to prevent path issues (e.g., slashes)
  sample=$(echo "$sample" | tr '/ ' '__')
  
  if [[ "$command" == "annotate" ]]; then
    pathToReads=$(pwd)/results/${sample}/${locus}/reads/${sample}_atleast-2.fasta
  fi

  # throttle: wait if we already have $max_jobs running
  while [ "$(squeue -u "$USER" -h | wc -l)" -ge "$max_jobs" ]; do
    sleep 1
  done
  
  # Create batch script content
  batch_script="#!/usr/bin/env bash
#SBATCH -p ${partition}
#SBATCH -J ${sample}_${command}
#SBATCH -o slog/${sample}_${command}.slog
#SBATCH --cpus-per-task=${cpus_per_task}
#SBATCH --exclusive

# skip Nextflow's internet/version check
export NXF_OFFLINE=1
export NXF_VER=22.10.6
export NXF_OPTS="-XX:ActiveProcessorCount=\$SLURM_CPUS_PER_TASK"
# module load nextflow   # if you need a module; otherwise remove

nextflow run ${NXF_SCRIPT} -offline \\
  -profile            \"$runtime\" \\
  --sample_name       \"$sample\" \\
  --reads             \"$pathToReads\" \\
  --outdir            \"./results/${sample}/${locus}\" \\
  --locus             $locus \\
  --species           Homo_sapiens \\
  --germline_ref_dir  \"$germline_ref_dir\" \\
  -with-report        \"./results/${sample}/${locus}/${sample}_nextflow_${locus}_${command}.html\" \\
  $additional_params"

  # If echo mode, print the commands
  if [[ "$echo_only" == true ]]; then
    echo "Would submit the following batch script for sample $sample:"
    echo "$batch_script"
    echo "---"
  else
    # the following cd addresses a weird java-related problem in WSL, which causes a failure with the error 
    # sbatch: error: getcwd failed: No such file or directory
    cd .
    # submit the batch script via heredoc
    sbatch <<< "$batch_script"
  fi

done < "$input_file"
