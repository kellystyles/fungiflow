#!/bin/bash
#SBATCH --array=78
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --partition=parallel
#SBATCH --time=4-00:00:00
#SBATCH --job-name=fungiflow
#SBATCH -o /nfs/scratch/styleske/fungiflow/projects/fungiflow/fungiflow_%A_%a.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kelly.styles@vuw.ac.nz

module load singularity/3.7.3

python3 fungiflow.py -d "test" -a ${SLURM_ARRAY_TASK_ID} \
-if /nfs/scratch/styleske/fungiflow/projects/fungiflow/test/${SLURM_ARRAY_TASK_ID}_F.fq.gz \
-ir /nfs/scratch/styleske/fungiflow/projects/fungiflow/test/${SLURM_ARRAY_TASK_ID}_R.fq.gz \
-n /nfs/scratch/styleske/fungiflow/projects/fungiflow/test/${SLURM_ARRAY_TASK_ID}_ONT_reads.fq.gz \
-c ${SLURM_CPUS_PER_TASK} -m ${SLURM_MEM_PER_NODE} -ant -t "hybrid" -dt "isolate" -b \
-s /nfs/scratch/styleske/fungiflow/projects/fungiflow/singularity_images/fungiflow_v0.2.sif \
-sf /nfs/scratch/styleske/fungiflow/projects/fungiflow/singularity_images/funannotate_gm.sifcp \
-sa /nfs/scratch/styleske/fungiflow/projects/fungiflow/singularity_images/antismash_6.0.1--pyhdfd78af_0.sif \
-idb /nfs/scratch/styleske/fungiflow/databases/fungi_ITS/ITS_RefSeq_Fungi -kdb /nfs/scratch/styleske/fungiflow/databases/kraken2_std/ \
-bdb /nfs/scratch/styleske/fungiflow/databases/NCBI_nt/nt
