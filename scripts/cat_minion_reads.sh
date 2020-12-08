#!/bin/bash
#SBATCH --cpus-per-task=12
#SBATCH --mem=24G
#SBATCH --partition=parallel
#SBATCH --time=1-00:00:00
#SBATCH --job-name=cat_minion_reads
#SBATCH -o /nfs/scratch/styleske/fungiflow/scripts/cat_minion_reads_%A_%a.out
#SBATCH -e /nfs/scratch/styleske/fungiflow/scripts/cat_minion_reads_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kelly.styles@vuw.ac.nz

# shouldn't need to trim Nanopore reads as the guppy software automatically does this when it calls the bases from the MinION output
(date && echo "| concatenating all nanopore reads...")
cat "$TWD"/data/raw/long/${SLURM_ARRAY_TASK_ID}_*.fastq > ${SLURM_ARRAY_TASK_ID}_minion.fq
