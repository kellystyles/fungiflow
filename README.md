# fungiflow

## Overview

Repeatable, parallelized, readable, and user friendly workflow for identifying fungal biosynthetic gene clusters (BGCs) from raw Illumina sequence data with minimal inputs. This pipeline is mostly written in `python3`.

The overall workflow is defined by several modules:
1. Assembly module
    1. Trimming of short Illumina reads with Trimmomatic
    2. Trimming of long MinION reads with porechop
    3. Correction of trimmed MinION reads with FMLRC, using trimmed short reads
    4. (OPTIONAL) Filtering out non-eukaroytic contamination with Kraken2 (standard database | Oct 2020)
    5. QC of trimmed and filtered reads with FastQC
    6. Assembly:
        - SPAdes for assembly with short Illumina sequence reads
        - Flye for hybrid assembly with short Illumina and long MinION sequence reads
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

![Overview of fungiflow pipeline](https://github.com/kellystyles/fungiflow/blob/main/workflow.png)
Overview of fungiflow pipeline

## Installation
Clone this GitHub repositroy by entering ```git clone https://github.com/kellystyles/fungiflow.git``` in the directory you would like to install in. You could add this directory to your PATH. 
```
export PATH=$PATH:/path/to/fungiflow/ >> ~/.bashrc
source ~/.bashrc
```

This pipeline runs from several Singularity containers to ensure a repeatable and cosistent output from input data. Consequently, the only dependency required is `singularity`. Some third party `python3` libraries are required, including `numpy` and `pandas`, and optionally `seaborn`. These can be easily installed using `conda` as so:
```
conda create -n fungiflow python3
source activate fungiflow
conda install numpy pandas seaborn
conda deactivate
```

This workflow is designed to operate on an HPC, so a lot of cpus and memory is required. I would suggest a minimum of 8 cpus and 40 GB of memory, as some of the programs in the workflow require a lot of memory.

If you wish to use the optional Funannotate module, you will need to obtain a new license to use GeneMark-ES software, obtainable from [here](http://topaz.gatech.edu/GeneMark/license_download.cgi). Extract the `gm_key_64` file to `~/` as so:
'''
tar -xvzf gm_key_64.tar.gz -O gm_key & mv gm_key ~/.gm_key
'''

~To check that all the necessary files are present the `check_install.sh` script can be executed from `/fungiflow/scripts`.~

### Databases

Several databases are needed for optional fungiflow modules to operate correctly. These can be installed by running the `install.py` script.
```
python3 install.py --help
```
All core databases are bundled within the Singularity containers so you can perform a basic run straight out of the box. The following optional databases are pretty large so the install may take several hours to download with 4 CPUs:
    - Kraken2 standard database ( 213 GB | ~6 hours); required for taxonomic filtering of short reads.
    - ITS_Refseq_Fungi database ( 162 Mb | 20 mins); required for BLASTn of extracted ITS sequences.
    - NCBI-nt database (132 GB | ~4 hours); required for assigning taxonomy to contigs for the Blobplots module.

## Usage

Usage can be shown using `--help`.
```
python3 fungiflow.py --help

usage: fungiflow.py [-h] -d DIRECTORY -if
                    ILLUMINA_F -ir ILLUMINA_R -a     
                    ARRAY -c CPUS -m MEM [-ant]      
                    [-f] [-s SINGULARITY_IMAGE]      
                    [-data DATABASE_PATH]
                    [-sf SINGULARITY_FUNANNOTATE]    
                    [-sa SINGULARITY_ANTISMASH]      
                    [-its] [-k] [-e] [-b]
                    [-idb ITS_DB] [-kdb KRAKEN2_DB]  
                    [-bdb BLOB_DB] [-edb EGGNOG_DB]  
                    [-n NANOPORE]
                    [-t {isolate,metagenomic}]       
                    [-minlen MINIMUM_LENGTH]
                    [-mtm MIN_TRAINING_MODELS]       
                    [--careful] [--dry_run DRY_RUN]  

Fungiflow - the automated eukaryotic genomic
pipeline for fungi.

optional arguments:
  -h, --help            show this help message and   
                        exit
  -d DIRECTORY, --directory DIRECTORY
                        Working directory path.      
  -if ILLUMINA_F, --illumina_f ILLUMINA_F
                        Illumina short forward       
                        reads path
  -ir ILLUMINA_R, --illumina_r ILLUMINA_R
                        Illumina short reverse       
                        reads path
  -a ARRAY, --array ARRAY
                        a unique value
  -c CPUS, --cpus CPUS  Number of threads to use     
  -m MEM, --mem MEM     Amount of memory to use (in  
                        GB)
  -ant, --antismash     Add this argument if you     
                        would like to search the     
                        assembly for BGCs
  -f, --funannotate     Add this argument if you     
                        would like to annotate the   
                        assembly
  -s SINGULARITY_IMAGE, --singularity_image SINGULARITY_IMAGE
                        Primary Singularity image    
                        for Fungiflow
  -data DATABASE_PATH, --database_path DATABASE_PATH 
                        Path to installed databases  
  -sf SINGULARITY_FUNANNOTATE, --singularity_funannotate SINGULARITY_FUNANNOTATE
                        Singularity image for        
                        Funannotate
  -sa SINGULARITY_ANTISMASH, --singularity_antismash SINGULARITY_ANTISMASH
                        Singularity image for        
                        antiSMASH
  -its, --its           Search assembly for ITS      
                        sequences. Part of post-     
                        analysis module.
  -k, --kraken2         Run trimmed reads against    
                        Kraken2 standard database.   
                        Will save unclassified       
                        reads (i.e., not matching    
                        the standard database),      
                        which will be used for       
                        assembly. Part of taxonomic  
                        module.
  -e, --eggnog          Functionally annotate the    
                        assembly proteins with       
                        eggnog. Part of annotation   
                        module.
  -b, --blobplot        Run blobtools module on      
                        output assembly. Part of     
                        blobtools module.
  -idb ITS_DB, --its_db ITS_DB
                        Path to alternative
                        ITS_refseq BLASTn database.  
  -kdb KRAKEN2_DB, --kraken2_db KRAKEN2_DB
                        Path to alternative Kraken2  
                        standard database.
  -bdb BLOB_DB, --blob_db BLOB_DB
                        Path to alternative NCBI-nt  
                        database for blobtools.      
  -edb EGGNOG_DB, --eggnog_db EGGNOG_DB
                        Path to alternative eggnog   
                        database for blobtools.      
  -n NANOPORE, --nanopore NANOPORE
                        Path to Nanopore MinION reads.
  -t {isolate,metagenomic}, --type {isolate,metagenomic}
                        Sequence data source type.   
                        Accepted arguments are       
                        'isolate' and 'metagenomic'  
  -minlen MINIMUM_LENGTH, --minimum_length MINIMUM_LENGTH
                        Minimum length of long-      
                        reads to retain during QC    
                        processes. Default is 2000   
                        bp.
  -mtm MIN_TRAINING_MODELS, --min_training_models MIN_TRAINING_MODELS
                        Mininum number of predicted  
                        models to train Augustus.    
                        Default is 200.
  --careful             Assembles reads with SPAdes  
                        using lower k-mer values
                        and runs in single cell      
                        mode.
  --dry_run DRY_RUN     Perform a dry run of the     
                        pipeline and check the       
                        input files, databases and   
                        expected output files 
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
│   │   (array1_val)_adapters.fasta       # FASTA file of adapters found in reads
│   │   (array1_val)_adapters.fasta
│   │   ...
│
└───trimmed
│   │   (array1_val)_trimmed_1P.fq.gz     # short paired forward trimmed reads
│   │   (array1_val)_trimmed_1U.fq.gz     # short unpaired forward trimmed reads
│   │   (array1_val)_trimmed_2P.fq.gz     # short paired reverse trimmed reads
│   │   (array1_val)_trimmed_2U.fq.gz     # short unpaired reverse trimmed reads
│   │   (array1_val)_ONT_reads.fq         # trimmed MinION reads
│   │   (array1_val)_ONT_reads.fa         # trimmed MinION reads in FASTA format
│   │   (array1_val)_ONT_reads_lf.fa      # length-filtered MinION reads
│   │   (array1_val)_ONT_reads_corr.fa    # corrected MinION reads
│   │   (array1_val)_mapped.npy           # short reads mapped to MinION reads
│   │   ...
│
└───kraken2
│   │   (array1_val)_class_1.fq           # classified forward reads
│   │   (array1_val)_unclass_1.fq         # unclassified forward reads
│   │   (array1_val)_class_2.fq           # classified reverse reads
│   │   (array1_val)_unclass_2.fq         # unclassified reverse reads
|   |   ...
│
└───assembly
│   │   (array1_val)_scaffolds.fasta      # SPADes assembly file
|   |   assembly.fasta                    # Flye assembly file
|   |   (array1_val)_ONT_corr.sam         # short reads mapped to Flye assembly
│   │   (array1_val)_pilon.fa             # pilon polished hybrid assembly file
│   │   racon_consensus.fa                # racon polished hybrid assembly file
|   |   racon_consensus.fa.bam            # short reads mapped to racon assembly
|   |   racon_consensus.fa.bam.bai        # index of above file
|   |   racon_consensus.hdf               # corrected MinION reads aligned to racon assembly
|   |   medaka_consensus.fa               # medaka polished hybrid assembly file
|   |   medaka_consensus.sam              # corrected MinION reads aligned to medaka assembly 
|   |   medaka_consensus.bai              # index of above file
|   |   medaka_sorted.sam                 # sorted file of `medaka_consensus.sam` 
|   |   medaka_sorted.sam.bai             # index of above file
|   |   polypolish.fasta                  # final polished hybrid assembly
│   
└───funannotate
│   │
│   └───(array1_val)
│       └───logfiles
│       └───predict_misc
│       └───predict_results
│       |   |   (array1_val).gbk          # GBK file of aseembly with predicted genes
|       └───eggnog
|       |   |   (array1_val).emapper.annotations  # eggnog annotations
│       └───annotate_misc
|       |   |   ...
│       └───annotate_results 
|           |   (array1_val).gbk          # GBK file with annotated genes       
│           |   ...
|
└───antismash
│   │
│   |   (array1_val).gbk                  # GBK file input to antiSMASH
|   |   (array1_val).json                 # JSON file of antiSMASH output
|   |   index.html                        # output HTML viewer
|   |   ...
|
└───quast
│   │
│   |   transposed_report.tsv             # output report file used by this pipeline
|   |   report.html                       # output HTML viewer
|   |   ...
|
└───blobplots
│   │
│   |   assembly.bam                      # short reads mapped to assembly
|   |   assembly.bam.bai                  # index of above file
|   |   assembly.blobDB.json              # blobtools JSON output file
|   |   megablast.out                     # MegaBLAST output
|   |   ...

```   
### Speed
Currently, fungiflow is designed for use in an HPC environment so will run quicker with more resources assigned to it. Below are some runtimes and hardware usage for a selection of fungal genomes.

## Planned implementations

- Work will be done to implement multiprocessing for slower parts of the pipeline, particularly lookup/identification tasks (e.g., `MegaBLAST` in the blobplots package).
- Implementation of assembly using only MinION long reads, particularly with the release of the R10 flow cells which purport a >99% accuracy rate.
- Whilst the repeatability and accessibility is ensured by the usage of Singularity containers, a conda environment and PyPI package is planned.