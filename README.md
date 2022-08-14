# fungiflow

Repeatable, parallelized, readable, and user friendly workflow for identifying biosynthetic gene clusters (BGCs) from raw short or long sequence data with minimal inputs - focused on fungal sequences.

The overall workflow is defined by several modules:
1. Assembly module
    1. Trimming of short Illumina reads with Trimmomatic
    2. Trimming of long ONT reads with porechop
    3. Correction of trimmed ONT reads with FMLRC, using trimmed short reads
    4. (OPTIONAL) Filtering out non-eukaroytic contamination with Kraken2 (standard database | Oct 2020)
    5. QC of trimmed and filtered reads with FastQC
    6. Assembly:
        - SPAdes for assembly with short Illumina sequence reads
        - Flye for hybrid assembly with short Illumina and long ONT sequence reads
2. (OPTIONAL) Annotation module:
    1. Cleaning, sorting, soft masking, and gene prediction of contigs with the Funannotate pipeline (v1.8.3)
3. Post-analysis module:
    1. Evaluations of assemblies using QUAST
    2. (OPTIONAL) Annotation and extraction of ITS regions using ITSx and BLASTn of ITS sequences (ITS_Refseq_Fungi | Oct 2020)
    3. (OPTIONAL) Identification of BGCs using antiSMASH (v5.1.3)
    4. Report generation
4. Blobplot module:
    1. Minimap2 maps short reads to assembly file in BAM format
    2. Samtools indexes the BAM file 
    3. MegaBLAST assigns taxonomy to each contig (NCBI-nt | Jun 2022)
    4. Blobtools generates blobplots

![Overview of fungiflow pipeline](https://github.com/kellystyles/fungiflow/blob/main/Pipeline_progress_Jan_2021.png)
Overview of fungiflow pipeline

## Installation
Clone this GitHub repositroy by entering ```git clone https://github.com/kellystyles/fungiflow.git``` in the directory you would like to install in.

This workflow is designed to operate on an HPC, so a lot of cpus and memory is required. I would suggest a minimum of 8 cpus and 40 GB of memory, as some of the programs in the workflow require a lot of memory.

To make this workflow repeatable, readable, and minimally difficult to operate, all the software needed has been packaged into a series of Singularity images. The only pre-requisite you need for the workflow to function other than a bash shell is Singularity.

- bash environment (tested with Ubuntu 18.02)
- Singularity >/ 3.5.3

If you wish to use the optional Funannotate module, you may need to copy the `gm_key` file to `~/.gm_key` on the cluster or machine you are working on. You may need to obtain a new license from [here](http://topaz.gatech.edu/GeneMark/license_download.cgi). Extract the `gm_key_64` file to `~/`

~To check that all the necessary files are present the `check_install.sh` script can be executed from `/fungiflow/scripts`.~

### Databases

Several databases are needed for optional fungiflow modules to operate correctly. One or all can be installed by running the `install.py` script.
```
python3 install.py --help
```
Some of these databases are pretty large so the install may take several hours:
    - Kraken2 standard database ( 213 GB | ~6 hours with 4 CPUs); required for taxonomic filtering of short reads.
    - ITS_Refseq_Fungi database ( 162 Mb | 20 mins with 4 CPUs); required for BLASTn of extracted ITS sequences.
    - NCBI-nt database (132 GB | ~4 hours with 4 CPUs); required for assigning taxonomy to contigs for the Blobplots module.

## Usage

Usage can be shown using `--help`.
```
python3 fungiflow.py --help

usage: fungiflow.py [-h] -d DIRECTORY -if ILLUMINA_F -ir ILLUMINA_R -a ARRAY
                    -c CPUS -m MEM [-ant] [-s SINGULARITY]
                    [-sf SINGULARITY_FUNANNOTATE] [-sa SINGULARITY_ANTISMASH]
                    [-idb ITS_DB] [-kdb KRAKEN2_DB] [-bdb BLOB_DB]
                    [-n NANOPORE] [-t {isolate,metagenomic}]
                    [-minlen MINIMUM_LENGTH] [-b]

Fungiflow - the automated eukaryotic genomic pipeline for fungi.

optional arguments:
  -h, --help            show this help message and exit
  -d DIRECTORY, --directory DIRECTORY
                        Working directory path.
  -if ILLUMINA_F, --illumina_f ILLUMINA_F
                        Illumina short forward reads path
  -ir ILLUMINA_R, --illumina_r ILLUMINA_R
                        Illumina short reverse reads path
  -a ARRAY, --array ARRAY
                        a unique value
  -c CPUS, --cpus CPUS  Number of threads to use
  -m MEM, --mem MEM     Amount of memory to use (in GB)
  -ant, --antismash     Add this argument if you would like to search the
                        assembly for BGCs
  -s SINGULARITY, --singularity SINGULARITY, --singularity SINGULARITY
                        Primary Singularity image for Fungiflow
  -sf SINGULARITY_FUNANNOTATE, --singularity_funannotate SINGULARITY_FUNANNOTATE
                        Singularity image for Funannotate
  -sa SINGULARITY_ANTISMASH, --singularity_antismash SINGULARITY_ANTISMASH
                        Singularity image for antiSMASH
  -idb ITS_DB, --its_db ITS_DB
                        Path to ITS_refseq BLASTn database. Part of taxonomic
                        module.
  -kdb KRAKEN2_DB, --kraken2_db KRAKEN2_DB
                        Path to Kraken2 standard database. Part of taxonomic
                        module.
  -bdb BLOB_DB, --blob_db BLOB_DB
                        Path to NCBI-nt database for blobtools. Part of
                        blobtools module.
  -n NANOPORE, --nanopore NANOPORE
                        Nanopore reads.
  -t {isolate,metagenomic}, --type {isolate,metagenomic}
                        Sequence data source type. Accepted arguments are
                        'isolate' and 'metagenomic'
  -minlen MINIMUM_LENGTH, --minimum_length MINIMUM_LENGTH
                        Minimum length of long-reads to retain during QC
                        processes. Default is 2000 bp.
  -b, --blobplot        Prepare blobplots of assembly
```

There is a SLURM script for running the pipeline on your SLURM-compatible HPC. Edit this with your specific variables prior to use.
``` 
sbatch fungiflow_slurm.sh
```

### Files/Directory Tree
After cloning into the fungiflow GitHub repository, create a new project folder in `/fungiflow/projects/`. In this folder, create a nested directory named `data/raw` and place all raw Illumina sequence reads into this folder in `*fq.gz format` (no preprocessed reads). The final path should look like `/fungiflow/projects/project_name/data/raw`.

To utilise parallel computing, sequence read names must be appeneded with an array value and pair information. This can be automatically performed by running `name_change.sh` in the `/fungiflow/scripts` folder as below. This will return a `filenames.txt` document which lists the name of each sequence file and the number prepended to it.
```
bash name_change.sh ...
```

```
project directory
│   (array1_val)_F.fq.gz
│   (array1_val)_R.fq.gz
│   (array1_val)_ONT_reads.fq.gz  
│   ...
│
└───adapters
│   │   (array1_val)_adapters.fasta
│   │   (array1_val)_adapters.fasta
│   │   ...
│
└───trimmed
│   │   (array1_val)_trimmed_1P.fq.gz     # short paired forward trimmed reads
│   │   (array1_val)_trimmed_1U.fq.gz     # short unpaired forward trimmed reads
│   │   (array1_val)_trimmed_2P.fq.gz     # short paired reverse trimmed reads
│   │   (array1_val)_trimmed_2U.fq.gz     # short unpaired reverse trimmed reads
│   │   (array1_val)_ONT_reads.fq         # trimmed ONT reads
│   │   (array1_val)_ONT_reads.fa         # trimmed ONT reads in FASTA format
│   │   (array1_val)_ONT_reads_lf.fa      # length-filtered ONT reads
│   │   (array1_val)_ONT_reads_corr.fa    # corrected ONT reads
│   │   (array1_val)_mapped.npy           # short reads mapped to ONT reads
│   │   ...
│
└───kraken2
│   │   (array1_val)_class_1.fq           # classified forward reads
│   │   (array1_val)_unclass_1.fq         # unclassified forward reads
│   │   (array1_val)_class_2.fq           # classified reverse reads
│   │   (array1_val)_unclass_2.fq         # unclassified reverse reads
│
└───assembly
│   │
│   └───(array1_val)
│       │   (array1_val)_scaffolds.fasta  # SPADes assembly file
│       │   (array1_val)_scaffolds.fasta  # Flye assembly file
│       │   (array1_val)_pilon.fasta      # pilon polished hybrid assembly file
│       │   
│       └───racon
│       │   |   (array1_val)_racon.fasta  # racon polished hybrid assembly file
│       │   
│       └───medaka
│           |   (array1_val)_medaka.fasta # medaka polished hybrid assembly file
│   
└───funannotate
│   │
│   └───(array1_val)
│       └───logfiles
│       └───predict_misc
│       └───predict_results
│           |   (array1_val).gbk          # GBK file of aseembly with predicted genes
│           |   ...
│
└───antismash
│   │
│   └───(array1_val)
└───quast
│   │
│   └───(array1_val)
└───blobplots
│   │
│   └───(array1_val)
```   