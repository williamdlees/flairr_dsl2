# Run this command from the directory above your 'results' directoy, when all samples have been processed
# setup paths to nextflow scripts
SCRIPT_PATH="${BASH_SOURCE[0]}"
# Convert to absolute path
SCRIPT_PATH="$(readlink -f "$SCRIPT_PATH")"
SCRIPT_DIR="$(dirname "$SCRIPT_PATH")"
NF_ROOT="${SCRIPT_DIR%/*}"

singularity exec -B $(pwd):/data -B $NF_ROOT:/nf_root ~/sifs/flairr_dsl2_latest.sif bash /nf_root/processed_samples/run_summaries.sh
