# fungiflow

Repeatable, readable, and minimalist workflow for identifying biosynthetic gene clusters from raw sequence data with minimal inputs with a focus on fungal sequences.

The overall workflow is defined by:
1. QC and trimming of reads with FastQC and Trimmomatic
2. Filtering non-eukaroytic contamination with Kraken2 (standard database | Oct 2020)
3. QC of trimmed and filtered reads with FastQC
4. Assembly:
    - SPAdes for single organisms sequence reads
    - metaSPAdes for metagenomic sequence reads
    - MaSurCa is used if the first assembly with SPAdes or metaSPAdes fails
5. Taxonomic identification of contigs using Kraken2 (database of all representative fungal genomes | Nov 2020)
6. Annotation an extraction of ITS regions using ITSx and BLASTn of ITS sequences (ITS_Refseq_Fungi | Oct 2020)
7. Evaluations of assemblies using QUAST
8. Identification of BGCs using antiSMASH (ver 5.1.3)

![Overview of fungiflow pipeline](https://github.com/kellystyles/fungiflow/blob/main/fungiflow_nov_2020.png)
Overview of fungiflow pipeline

## Installation
Clone this GitHub repositroy by entering ```git clone https://github.com/kellystyles/fungiflow.git``` in the directory you would like to install in.

This workflow is designed to operate on an HPC, so a lot of cpus and memory is required. I would suggest a minimum of 8 cpus and 40 GB of memory, as some of the programs in the workflow require a lot of memory.

To make this workflow repeatable, readable, and minimally difficult to operate, all the software needed has been packaged into a series of Singularity images. The only pre-requisite you need for the workflow to function other than a bash environment is Singularity.

- bash environment (tested with Ubuntu 18.02)
- Singularity >/ 3.5.3

## Usage

After cloning into the fungiflow GitHub repository, create a new project folder in `/fungiflow/projects/`. In this folder, create a nested directory named `data/raw` and place all raw Illumina sequence reads into this folder. The final path should look like `/fungiflow/projects/project_name/data/raw`.

To utilise parallel computing, the sequence read names must be appeneded with a number and pair information. This can be automatically performed by running `name_change.sh` in 
the `/fungiflow/scripts` folder as below.
```
name_change.sh /fungiflow/projects/project_name/data/raw
```

Following this, ensure that all the Singularity images are present in the `/fungiflow/images/` folder:
1. qc.sif
2. assembly.sif
3. kraken2.sif
4. itsx.sif
5. antismash.sif

Make a copy of fungiflow.sh
```
cp fungiflow.sh fungiflow_copy.sh
```
Edit this copy with the desired arguments under `### Variables ###`, as well as in the SLURM `#SBATCH` arguments. These arguments must match (e.g., if `#SBATCH --cpus-per-task=12` then the following is required `cpus=12`). Ensure that the correct number of tasks has been designated in the SLURM array (e.g., `#SBATCH --array=1-n%x`; where `n` is the total number of tasks and `x` is the number of tasks to run at once).

Submit this script to SLURM
``` 
sbatch fungiflow_copy.sh
```

The resulting report can be found at `/fungiflow/projects/project_name` directory for each individual sequence input and will be labelled with the format: `<time>_<project_name>_fungiflow_<task_id>.txt`.
The output file can be found in the `/fungiflow/scripts` with the `<job_number>_<task_id>_<job_name>.log` naming convention.
