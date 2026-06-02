#!/bin/bash
#SBATCH -c 1                              # Request cores
#SBATCH -t 1-00:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=4G                       # Memory total in MiB (for all cores)
#SBATCH -o hostname_%j.out                 # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e hostname_%j.err                 # File to which STDERR will be written, including job ID (%j)
                                           # You can change the filenames given with -o and -e to any filenames you'd like
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=praefrontalis@gmail.com
########################################################################
# Argument parser
########################################################################
OPTIND=1         # Reset in case getopts has been used previously in the shell.
# Initialize our own variables:
WORK_DIR="/n/scratch/users/w/wit498/workdir"
SAMPLE_NAME=""
NF_PARAMS=""
while getopts "w:s:" opt; do
  case "$opt" in
    w)  WORK_DIR=$OPTARG
      ;;
    s)  SAMPLE_NAME=$OPTARG
      ;;
  esac
done
shift $((OPTIND-1))
[ "${1:-}" = "--" ] && shift
NF_PARAMS=$@
echo "WORK_DIR=$WORK_DIR, SAMPLE_NAME=$SAMPLE_NAME, NF_PARAMS:$NF_PARAMS"

module load java/jdk-23.0.1 conda/miniforge3/24.11.3-0

mamba activate rnaseq

export NXF_SINGULARITY_CACHEDIR=/n/app/containers/shared/nf-core/rnaseq/3.14.0

nextflow run \
    nf-core/rnaseq \
    -r 3.14.0 \
    -c /home/wit498/projs/bedRest/nextflow-slurm.config \
    -profile singularity,cluster \
    -params-file /home/wit498/projs/bedRest/params.yaml \
    -w $WORK_DIR \
    -resume