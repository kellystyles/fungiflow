#!/bin/bash
#SBATCH --array=1-4%4
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --partition=parallel
#SBATCH --time=8:00:00
#SBATCH --job-name=fungiflow
#SBATCH -o /output/path/fungiflow_%A_%a.log
#SBATCH --mail-type=ALL


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
    If you wanted to use a specific name instead of a number, remove the 
    `#SBATCH --array` shebang and uncomment the `SLURM_ARRAY_TASK_ID` and 
    edit the text to the desired name. The input filename structure will still
    need to follow the same naming convention as for an array run.
    This can also be changed once everything is completed so its not biggie
 - `-t` runs pipeline in `isolate` mode
"""

# Uncomment if you want to use a name instead of a SLURM array number
#SLURM_ARRAY_TASK_ID=name

# load a conda and Singularity module
module load Miniconda3
module load singularity

# activate the fungiflow conda environment
source activate fungiflow

# perform a minimal run (hybrid assembly only)
# add other parameters as you like
python3 /path/to/fungiflow/fungiflow.py \
-d "path/to/directory/with/reads" \
-a ${SLURM_ARRAY_TASK_ID} \
-if ${SLURM_ARRAY_TASK_ID}_1.fq.gz \
-ir ${SLURM_ARRAY_TASK_ID}_2.fq.gz \
--nanopore ${SLURM_ARRAY_TASK_ID}_ONT_reads.fastq \
-c ${SLURM_CPUS_PER_TASK} -m ${SLURM_MEM_PER_NODE} -t "isolate" \
-s /path/to/fungiflow_v3.sif