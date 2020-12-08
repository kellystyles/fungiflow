#!/bin/bash
#SBATCH --mem=12G
#SBATCH --constraint=Intel
#SBATCH --partition=parallel
#SBATCH --time=12:00:00
#SBATCH --job-name=anti_concat_json
#SBATCH -o /nfs/scratch/styleske/fungiflow/scripts/anti_concat_json_%A_%a.out

module load singularity/3.5.2

singularity exec $1 \
run_antismash.py --reuse-results $2 --output-dir $3
