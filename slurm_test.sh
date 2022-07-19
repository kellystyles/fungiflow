#!/bin/bash
#SBATCH --array=78
#SBATCH --cpus-per-task=24
#SBATCH --mem=150G
#SBATCH --partition=parallel
#SBATCH --time=2-00:00:00
#SBATCH --job-name=fungiflow
#SBATCH -o /nfs/scratch/styleske/fungiflow/projects/fungiflow/fungiflow_%A_%a.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kelly.styles@vuw.ac.nz

module load singularity/3.7.3

python3 ./fungiflow.py -d "/nfs/scratch/styleske/fungiflow/projects/microcera_test" -a ${SLURM_ARRAY_TASK_ID} \
-if ${SLURM_ARRAY_TASK_ID}_F.fq.gz \
-ir ${SLURM_ARRAY_TASK_ID}_R.fq.gz \
-n ${SLURM_ARRAY_TASK_ID}_ONT_reads.fq.gz \
-c ${SLURM_CPUS_PER_TASK} -m ${SLURM_MEM_PER_NODE} -ant -t "isolate" -b \
-s /nfs/scratch/styleske/fungiflow/projects/singularity_files/fungiflow_v0-3_latest.sif \
-sf /nfs/scratch/styleske/fungiflow/projects/singularity_files/funannotate_gm.sifcp \
-sa /nfs/scratch/styleske/fungiflow/projects/singularity_files/antismash_6.0.1--pyhdfd78af_0.sif \
-idb /nfs/scratch/styleske/fungiflow/databases/fungi_ITS/ITS_RefSeq_Fungi -kdb /nfs/scratch/styleske/fungiflow/databases/kraken2_std/ \
-bdb /nfs/scratch/styleske/fungiflow/databases/NCBI_nt/nt
