#!/bin/bash
#
#SBATCH -J safe
#SBATCH -p short
#SBATCH -c 1
#SBATCH --time=01:00:00
#SBATCH --mem=2048M
#SBATCH -o logs/safe_%A.out
#SBATCH -e logs/safe_%A.err
#
cd /wrk/users/gyntartu/pipeline/out/uniprot_sprot.40.multi-step
n=$SLURM_ARRAY_TASK_ID
srun /proj/gyntartu/alignment-safety/safety-windows/main -f clean/cluster_${n}.clean.fasta > safety/cluster_${n}.out