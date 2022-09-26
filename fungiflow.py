#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import datetime
import library as lib
import assembly, blobplot, funannotate, post_analysis

"""
Main script of the Fungiflow pipeline.
USAGE: python3 fungiflow -d DIRECTORY -if ILLUMINA_F -ir ILLUMINA_R -a ARRAY -c CPUS -m MEMORY [EXTRA ARGUMENTS]
"""

def get_args():
    """Parse command line arguments"""
    
    try:
        parser = argparse.ArgumentParser(
            description="Fungiflow - the automated eukaryotic genomic pipeline for fungi.")
        parser.add_argument('-d', '--directory', action='store',
                            help='Working directory path.', type=str, required=True)     
        parser.add_argument('-if', '--illumina_f', action='store',
                            help='Illumina short forward reads path', type=str, required=True)
        parser.add_argument('-ir', '--illumina_r', action='store',
                            help='Illumina short reverse reads path', type=str, required=True)
        parser.add_argument('-a', '--array', action='store',
                            help='a unique value', type=str, required=True)
        parser.add_argument('-c', '--cpus', action='store',
                            help='Number of threads to use', type=str, required=True)
        parser.add_argument('-m', '--mem', action='store',
                            help='Amount of memory to use (in GB)', type=str, required=True)
        parser.add_argument('-ant', '--antismash', action='store_true',
                            help='Add this argument if you would like to search the assembly for BGCs')
        parser.add_argument('-f', '--funannotate', action='store_true',
                            help='Add this argument if you would like to annotate the assembly')
        parser.add_argument('-s', '--singularity_image', action='store',
                            help='Primary Singularity image for Fungiflow', type=str)
        parser.add_argument('-sf', '--singularity_funannotate', action='store',
                            help='Singularity image for Funannotate', type=str)
        parser.add_argument('-sa', '--singularity_antismash', action='store',
                            help='Singularity image for antiSMASH', type=str)
        parser.add_argument('-idb', '--its_db', action='store',
                            help='Path to ITS_refseq BLASTn database. Part of taxonomic module.', type=str)
        parser.add_argument('-kdb', '--kraken2_db', action='store',
                            help='Path to Kraken2 standard database. Part of taxonomic module.', type=str)
        parser.add_argument('-bdb', '--blob_db', action='store',
                            help='Path to NCBI-nt database for blobtools. Part of blobtools module.', type=str)
        parser.add_argument('-n', '--nanopore', action='store',
                            help='Nanopore reads.', type=str)
        parser.add_argument('-t', '--type', action='store', default="isolate", 
                            help='Sequence data source type. Accepted arguments are \
                            \'isolate\' and \'metagenomic\'', type=str, choices=['isolate', 'metagenomic'])
        parser.add_argument('-minlen', '--minimum_length', action='store', default="2000", 
                            help='Minimum length of long-reads to retain during QC processes. Default is 2000 bp.', type=str)                            
        parser.add_argument('-b', '--blobplot', action='store_true', 
                            help='Prepare blobplots of assembly')            
    except argparse.ArgumentError:
        lib.print_e("An exception occurred with argument parsing. Check your inputs.")
        parser.print_help(sys.stderr)
        sys.exit()

    return parser.parse_args()

### Begin Main Script ###
def main():
    input_args = get_args()
    lib.print_t("⁂⁂⁂⁂⁂⁂⁂⁂ Script Begins ⁂⁂⁂⁂⁂⁂⁂⁂  \n")
    start_time = datetime.datetime.now()
    print(input_args)

    ### MAKE INPUT DIR COMPLETE PATH
    input_args.directory_new = os.path.abspath(os.path.join(input_args.directory, input_args.array))
    lib.make_path(input_args.directory_new)
    os.chdir(input_args.directory_new)
    
    # Sets up class object for filenames
    filenames = lib.Files(input_args)

    # Running 'ASSEMBLY' module
    assembly.main(input_args, filenames)

    # Running 'FUNANNOTATE' module
    if input_args.funannotate is not None:
        funannotate.main(input_args, filenames)
    else:
        lib.print_h("Skipping assembly annotation...")
        filenames.funannotate_gbk = filenames.assembly_fasta

    # Running 'POST_ANALYSIS' module
    post_analysis.main(input_args, filenames)

    # Running 'BLOBPLOT' module
    if input_args.blobplot is not None:
        blobplot.main(input_args, filenames)

    lib.print_h(f"Script completed assembly and analysis of the sequence data in {datetime.datetime.now() - start_time}")
    lib.print_h(f"Results saved to {filenames.csv_output}")
    lib.print_h("Output files and variables are listed below:")
    filenames.printer()
    lib.print_tu("\n⁂⁂⁂⁂⁂⁂⁂⁂ Script Finished ⁂⁂⁂⁂⁂⁂⁂⁂")

if __name__ == '__main__':
    main()
