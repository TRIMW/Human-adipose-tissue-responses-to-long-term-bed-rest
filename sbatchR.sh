#!/bin/bash
#SBATCH -c 1                              # Request cores
#SBATCH -t 1-00:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=64G                       # Memory total in MiB (for all cores)
#SBATCH -o hostname_%j.out                 # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e hostname_%j.err                 # File to which STDERR will be written, including job ID (%j)
                                           # You can change the filenames given with -o and -e to any filenames you'd like
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=kseniia_petrova@hms.harvard.edu
########################################################################

module load gcc/14.2.0 R/4.4.2

Rscript /home/wit498/projs/Human-adipose-tissue-responses-to-long-term-bed-rest/1_differential_coexpression_for_cluster.R
