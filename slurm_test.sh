#!/bin/bash
#SBATCH --array=1-4%4
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --partition=parallel
#SBATCH --time=8:00:00
#SBATCH --job-name=fungiflow
#SBATCH -o /nfs/scratch/campbel2/some_directory_name/fungiflow_%A_%a.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kelly.styles@vuw.ac.nz


"""
This script will run the Fungiflow pipeline

Setup:
 - short reads
    must be in fq.gz format, with the suffix `_1.fq.gz` for forward reads
    and `_2.fq.gz` for reverse reads
 - ONT reads
    must be uncompressed and in a single fastq file
 - reads must be prefixed with SLURM_ARRAY_TASK_ID value if running in parallel
    i.e., if using `#SBATCH --array`

Paramterers:
 - `-a` is the array value that will be passed to the pipeline, and the files 
    will use this as a prefix. Here I have used the SLURM_ARRAY_TASK_ID as the 
    array value as then its easier to run multiple datasets in parallel. Both
    short reads and ONT reads should then be prefixed with the array value, e.g.
        1_1.fq.gz 1_2.fq.gz 1_ONT_reads.fastq 
        2_1.fq.gz 2_2.fq.gz 2_ONT_reads.fastq 
        ...
    If you wanted to use a specific name instead of a number I can show you how. 
    This can also be changed once everything is completed so its not biggie
 - `-t` runs pipeline in isolate mode
"""


module load old-mod-system/Miniconda3/4.9.2
module load singularity/3.7.3

source activate fungiflow

python3 /path/to/fungiflow/fungiflow.py \
-d "path/to/directory/with/reads" \
-a ${SLURM_ARRAY_TASK_ID} \
-if ${SLURM_ARRAY_TASK_ID}_1.fq.gz \
-ir ${SLURM_ARRAY_TASK_ID}_2.fq.gz \
--nanopore ${SLURM_ARRAY_TASK_ID}_ONT_reads.fastq \
-c ${SLURM_CPUS_PER_TASK} -m ${SLURM_MEM_PER_NODE} -t "isolate" \
-s /path/to/fungiflow_v3.sif