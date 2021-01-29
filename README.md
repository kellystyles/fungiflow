# fungiflow

Repeatable, parallelized, readable, and minimalist workflow for identifying biosynthetic gene clusters (BGCs) from raw sequence data with minimal inputs - focused on fungal genomic or metagenomic sequences.

The overall workflow is defined by:
1. QC and trimming of reads with FastQC and Trimmomatic
2. Filtering non-eukaroytic contamination with Kraken2 (standard database | Oct 2020)
3. QC of trimmed and filtered reads with FastQC
4. Assembly:
    - SPAdes for genomic sequence reads
    - metaSPAdes for metagenomic sequence reads
    - MaSurCa is used if the first assembly with SPAdes or metaSPAdes fails
5. Taxonomic identification of contigs using Kraken2 (database of all representative fungal genomes | Nov 2020)
6. Annotation an extraction of ITS regions using ITSx and BLASTn of ITS sequences (ITS_Refseq_Fungi | Oct 2020)
7. Evaluations of assemblies using QUAST
8. Identification of BGCs using antiSMASH (v5.1.3)
9. Annotation of contigs with the Funannotate pipeline (v1.8.3)

![Overview of fungiflow pipeline](https://github.com/kellystyles/fungiflow/blob/main/fungiflow_nov_2020.png)
Overview of fungiflow pipeline

## Installation
Clone this GitHub repositroy by entering ```git clone https://github.com/kellystyles/fungiflow.git``` in the directory you would like to install in.

This workflow is designed to operate on an HPC, so a lot of cpus and memory is required. I would suggest a minimum of 8 cpus and 40 GB of memory, as some of the programs in the workflow require a lot of memory.

To make this workflow repeatable, readable, and minimally difficult to operate, all the software needed has been packaged into a series of Singularity images. The only pre-requisite you need for the workflow to function other than a bash environment is Singularity.

- bash environment (tested with Ubuntu 18.02)
- Singularity >/ 3.5.3

If you wish to use the optional Funannotate module, you will need to copy the `gm_key` file to `~/.gm_key` on the cluster or machine you are working on. You could use a command like `cp /fungiflow/scripts/gm_key ~/.gm_key`. This key has an academic license and will only last for 400 days before it needs to be re-used. If this is the case, the GeneMark-ES software will fail (can check in the log file) and a new license must be downloaded from [here](http://topaz.gatech.edu/GeneMark/license_download.cgi). Extract the `gm_key_64` file to `~/`

## Usage

Most modules in the FungiFlow pipeline are optional and can be performed individually by calling the specific module in the parameters like so...
```singularity exec /path/to/fungiflow/images/module.sif module_name.sh --input file.fasta --output file_modulized --cpus 24```

### Files/Directory Tree
After cloning into the fungiflow GitHub repository, create a new project folder in `/fungiflow/projects/`. In this folder, create a nested directory named `data/raw` and place all raw Illumina sequence reads into this folder in `*fq.gz format` (no preprocessed reads). The final path should look like `/fungiflow/projects/project_name/data/raw`.

To utilise parallel computing, sequence read names must be appeneded with a number and pair information. This can be automatically performed by running `name_change.sh` in 
the `/fungiflow/scripts` folder as below. This will return a `filenames.txt` document which lists the name of each sequence file and the number prepended to it.
```
name_change.sh /fungiflow/projects/project_name/data/raw
```

### Databases

Several databases are needed for fungiflow to operate correctly. These can be installed by running the `install_db.sh` script as such from the `/fungiflow/scripts` directory
```
bash install_db.sh
```
Some of these databases are pretty large so the install may take several hours.

### Singularity images

Following this, ensure that all the Singularity images are present in the `/fungiflow/images/` folder:
1. qc.sif
2. assembly.sif
3. kraken2.sif
4. itsx.sif
5. antismash.sif
6. funannotate.sif

### Running fungiflow

Make a copy of fungiflow.sh
```
cp fungiflow.sh fungiflow_copy.sh
```
Edit this copy with the desired arguments under `### Variables ###`, as well as in the SLURM `#SBATCH` arguments. These arguments must match (e.g., if `#SBATCH --cpus-per-task=12` then the following is required `cpus=12`). Ensure that the correct number of tasks has been designated in the SLURM array (e.g., `#SBATCH --array=1-n%x`; where `n` is the total number of tasks and `x` is the number of tasks to run in parallel).

Submit this script to SLURM
``` 
sbatch fungiflow_copy.sh
```

## Output

The resulting report can be found at `/fungiflow/projects/project_name` directory for each individual sequence input and will be labelled with the format: `<TIME>_<PROJECT_NAME>_fungiflow_<SLURM_ARRAY_TASK_ID>.txt`. The output from `seff ${SLURM_JOB_ID}-${SLURM_ARRAY_TASK_ID}` will be appended to the end of this report to see resource utilisation.

 - the log file for each array task can be found in `/fungiflow/scripts` with the `<JOB_NUMBER>_<SLURM_ARRAY_TASK_ID>_<JOB_NAME>.log` naming convention.
 - adapters can be found in `/fungiflow/projects/<PROJECT_NAME>/data/adapters/` in 'fasta' format
 - trimmed reads can be found in `/fungiflow/projects/<PROJECT_NAME>/data/trimmed/` in 'fg.gz' format
 - the FastQC output files can be found in `/fungiflow/projects/<PROJECT_NAME>/output/(raw_QC|trimmed_QC)/`
 - Kraken2-filtered reads can be found in `/fungiflow/projects/<PROJECT_NAME>/data/classified/`
       - the classified reads are named `<SLURM_ARRAY_TASK_ID>_cseqs_1.fq` and `<SLURM_ARRAY_TASK_ID>_cseqs_2.fq`
       - the unclassified reads are named `<SLURM_ARRAY_TASK_ID>_useqs_1.fq` and `<SLURM_ARRAY_TASK_ID>_useqs_2.fq` (used for assembly)
 - the assembled reads can be found in `/fungiflow/projects/<PROJECT_NAME>/data/assembly/<SLURM_ARRAY_TASK_ID>/`
 - Kraken2 output of reads and contigs can be found in `/fungiflow/projects/<PROJECT_NAME>/output/classified/`
       - the reads output is named `<SLURM_ARRAY_TASK_ID>_kraken2_reads_report.txt`
       - the contig output is named `<SLURM_ARRAY_TASK_ID>_kraken2_contigs_report.txt`
       - the Kronogram for the contigs report is named `<SLURM_ARRAY_TASK_ID>_contigs_kronogram.html`
 - Quast output can be found in `/fungiflow/projects/<PROJECT_NAME>/output/quast/<SLURM_ARRAY_TASK_ID>/`
 - antiSMASH output can be found in `/fungiflow/projects/<PROJECT_NAME>/output/antismash/<SLURM_ARRAY_TASK_ID>/`
 - Funannotate output can be found in `/fungiflow/projects/<PROJECT_NAME>/output/funannotate/<SLURM_ARRAY_TASK_ID>/`
 
### Collating output and reports

For convenience, a script that will collate all the QC and output is available
```
bash report_collator.sh
```
This will create the following new files:
- report file at `/fungiflow/projects/<PROJECT_NAME>/output/<TIME>_<PROJECT_NAME>_report_collator.txt`
- MultiQC output in `/fungiflow/projects/<PROJECT_NAME>/output/<PROJECT_NAME>_(raw|trimmed)_multiqc_report.html`
- collated QUAST output in `/fungiflow/projects/<PROJECT_NAME>/output/collated_quast_reports.tsv`
- collated antiSMASH output in the `/fungiflow/projects/<PROJECT_NAME>/output/antismash_collated/` folder
- collate all annotations from the final Funannotate annotated genbank file into a tsv file in each annotation folder e.g., `/fungiflow/projects/<PROJECT_NAME>/output/funannotate/<SLURM_ARRAY_TASK_ID>/<SLURM_ARRAY_TASK_ID>_funannotate_annotations.tsv`
