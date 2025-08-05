#!/bin/bash
# process_slurm.sh
# Usage: ./process_slurm.sh <command> <input_fofn> <locus> <max_jobs> <container_runtime> [additional_nextflow_parameters...]
# <command> must be preprocess or annotate
# <input_fofn> should consist of lines of the form <sample_name>\t<fastq_path_name>
# <container_runtime> must be docker or singularity
# [additional_nextflow_parameters...] will be passed directly to nextflow
# slurm log files are written to the subdirectory slog

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

command="${1:?Usage: $0 <command> <input_tsv> <locus> <max_jobs> <container_runtime> [additional_nextflow_parameters...]}"
input_file="${2:?Usage: $0 <command> <input_tsv> <locus> <max_jobs> <container_runtime> [additional_nextflow_parameters...]}"
locus="${3:?Usage: $0 <command> <input_tsv> <locus> <max_jobs> <container_runtime> [additional_nextflow_parameters...]}"
max_jobs="${4:?Usage: $0 <command> <input_tsv>  <locus> <max_jobs> <container_runtime> [additional_nextflow_parameters...]}"
runtime="${5:?Usage: $0 <command> <input_tsv>  <locus> <max_jobs> <container_runtime> [additional_nextflow_parameters...]}"

# Shift the first 5 mandatory parameters
shift 5

# Any remaining parameters will be passed to nextflow
additional_params="$@"

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
    echo "$runtime is not available"
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

  

  # submit a little batch script via heredoc to avoid any nested-quoting issues
#sbatch <<EOF
#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -J ${sample}_${command}
#SBATCH -o slog/${sample}_${command}.slog
#SBATCH --cpus-per-task=12
#SBATCH --oversubscribe

# skip Nextflow's internet/version check
export NXF_OFFLINE=1
module load nextflow   # if you need a module; otherwise remove

echo nextflow run ${NXF_SCRIPT} -offline \
  -profile            "$runtime" \
  --sample_name       "$sample" \
  --reads             "$pathToReads" \
  --outdir            "./results/${sample}/IGH" \
  --locus             $locus \
  --species           Homo_sapiens \
  --germline_ref_dir  "/home/zmvanw01/ogr-ref" \
  $additional_params
#EOF

done < "$input_file"