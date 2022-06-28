#!/usr/bin/env
# -*- coding: utf-8 -*-

import os
import sys, argparse
import datetime
import library as lib
import assembly,blobplot,funannotate,post_analysis

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
        parser.add_argument('-s', '--singularity', '--singularity', action='store',
                            help='Primary Singularity image for Fungiflow', type=str, default=os.path.join("images","fungiflow.sif"))
        parser.add_argument('-sf', '--singularity_funannotate', action='store',
                            help='Singularity image for Funannotate', type=str, default=os.path.join("images","funannotate.sif"))
        parser.add_argument('-sa', '--singularity_antismash', action='store',
                            help='Singularity image for antiSMASH', type=str, default=os.path.join("images","antismash.sif"))
        parser.add_argument('-idb', '--its_db', action='store', default=os.path.join("databases","fungi_ITS_refseq","ITS"),
                            help='Path to ITS_refseq BLASTn database. Part of taxonomic module.', type=str)
        parser.add_argument('-kdb', '--kraken2_db', action='store', default=os.path.join("databases","kraken2_std"),
                            help='Path to Kraken2 standard database. Part of taxonomic module.', type=str)
        parser.add_argument('-bdb', '--blob_db', action='store', default=os.path.join("databases","NCBI_nt","nt"),
                            help='Path to NCBI-nt database for blobtools. Part of blobtools module.', type=str)
        parser.add_argument('-n', '--nanopore', action='store',
                            help='Nanopore reads.', type=str)
        parser.add_argument('-t', '--type', action='store', default="isolate", 
                            help='Sequence data source type. Accepted arguments are \
                            \'isolate\' and \'metagenomic\'', type=str, choices=['isolate', 'metagenomic'])
        parser.add_argument('-b', '--blobplot', action='store_true', 
                            help='Prepare blobplots of assembly')        

        if len(sys.argv) < 6:
            parser.print_help(sys.stderr)
            lib.print_e("Parser expecting 6 or more arguments")
            sys.exit()
    except argparse.ArgumentError:
        parser.print_help(sys.stderr)
        lib.print_e("An exception occurred with argument parsing. Check your inputs.")
        sys.exit()

    for arg in parser:
        print(f"{arg} \n")

    return parser.parse_args()

### Begin Main Script ###
def main():
    
    lib.print_tu("\n⁂⁂⁂⁂⁂⁂⁂⁂ Script Begins ⁂⁂⁂⁂⁂⁂⁂⁂\n")
    start_time = datetime.datetime.now()
    args = get_args()
    os.chdir(args.directory)

    # Sets up class object for filenames
    names = lib.Files(args.illumina_f,args.illumina_r)

    # Running 'ASSEMBLY' module
    assembly.main(args,names)

    # Running 'FUNANNOTATE' module
    if args.singularity_funannotate is not None:
        funannotate.main(args,names)
    else:
        lib.print_h("Skipping assembly annotation...")
        names.funannotate_gbk = names.assembly_fasta

    # Funannotate module changes into its own directory to prevent overlap of temp files with other array jobs
    # Checks the directory has changed back
    if os.getcwd() is not args.directory:
        lib.print_n("Changing back to main directory")
        os.chdir(args.directory)

    # Running 'POST_ANALYSIS' module
    post_analysis(args,names)

    # Running 'BLOBPLOT' module
    if args.blobplot is not None:
        blobplot.main(args,names)

    lib.print_h(f"Script completed assembly and analysis of the sequence data in {datetime.datetime.now() - start_time}")
    lib.print_h(f"Results saved to {names.csv_output}")
    lib.print_h("Output files are listed below:")
    lib.print_n(names)
    lib.print_tu("\n⁂⁂⁂⁂⁂⁂⁂⁂ Script Finished ⁂⁂⁂⁂⁂⁂⁂⁂")

if __name__ == '__main__':
    main()