#!/usr/bin/env python3
import os
import sys
import argparse
import datetime
import library as lib
import assembly, funannotate, post_analysis

"""
Main script of the Fungiflow pipeline.
USAGE: python3 fungiflow -d DIRECTORY -if ILLUMINA_F -ir ILLUMINA_R -a ARRAY -c CPUS -m MEMORY [OPTIONAL ARGUMENTS]
"""

def get_args():
    """Parse command line arguments"""
    
    try:
        parser = argparse.ArgumentParser(
            description="Fungiflow - the reproducible genomics pipeline for fungi.")
        parser.add_argument('-d', '--directory', action='store',
                            help='Working directory path.', type=str, required=True)     
        parser.add_argument('-a', '--array', action='store',
                            help='an informative and unique string', type=str, required=True)
        parser.add_argument('-if', '--illumina_f', action='store', default=None,
                            help='Illumina short forward reads path', type=str, required=False)
        parser.add_argument('-ir', '--illumina_r', action='store', default=None,
                            help='Illumina short reverse reads path', type=str, required=False)
        parser.add_argument('-n', '--nanopore', action='store', default=None,
                            help='Path to MinION reads.', type=str, required=False)
        parser.add_argument('-c', '--cpus', action='store',
                            help='Number of threads to use', type=str, required=True)
        parser.add_argument('-m', '--mem', action='store',
                            help='Amount of memory to use (in GB)', type=str, required=True)

        # Create mutually exclusive groups for various modules and their database paths
        # e.g., '-edb' cannot be supplied without '-e'
        parser.add_argument('-data', '--database_path', action='store', required=False, default=None,
                            help='Path to installed databases', type=str)
        parser.add_argument('-s', '--singularity_image', action='store', required=False,
                            help='Primary Singularity image for Fungiflow', type=str)
        parser.add_argument('-f', '--funannotate', action='store_true', required=False,
                            help='Add this argument if you would like to annotate the assembly')
        parser.add_argument('-sf', '--singularity_funannotate', action='store',
                            help='Singularity image for Funannotate', type=str, required=False)
        parser.add_argument('-ant', '--antismash', action='store_true', required=False,
                            help='Add this argument if you would like to search the assembly for BGCs')
        parser.add_argument('-sa', '--singularity_antismash', action='store',
                            help='Singularity image for antiSMASH', type=str, required=False)
        
        parser.add_argument('-e', '--eggnog', action='store_true', required=False,
                        help='Functionally annotate the assembly proteins with eggnog. \
                                Part of annotation module.')
        parser.add_argument('-edb', '--eggnog_db', action='store', required=False, 
                            help='Path to alternative eggnog database for eggnog.', type=str)
        parser.add_argument('-its', '--its', action='store_true', required=False,
                            help='Search assembly for ITS sequences. Part of post-analysis module.')
        parser.add_argument('-idb', '--its_db', action='store', required=False,
                            help='Path to alternative ITS_refseq BLASTn database.', type=str)
        parser.add_argument('-k', '--kraken2', action='store_true', required=False,
                            help='Run trimmed reads against Kraken2 standard database. \
                                Will save unclassified reads (i.e., not matching the standard database), \
                                    which will be used for assembly. Part of taxonomic module.')
        parser.add_argument('-kdb', '--kraken2_db', action='store', required=False,
                            help='Path to alternative Kraken2 standard database.', type=str)
        
        # optional arguments
        parser.add_argument('-t', '--type', action='store', default="isolate", required=False, 
                            help='Sequence data source type. Accepted arguments are \
                            \'isolate\' and \'metagenomic\'', type=str, choices=['isolate', 'metagenomic'])
        parser.add_argument('-minlen', '--minimum_length', action='store', default="2000", 
                            help='Minimum length of long-reads to retain during QC processes. \
                                Default is 2000 bp.', type=str, required=False)
        parser.add_argument('-mtm', '--min_training_models', action='store', default="200", 
                            help='Mininum number of predicted models to train Augustus. \
                                Default is 200.', type=str, required=False)
        parser.add_argument('--careful', action='store_true', required=False,
                            help='Assembles reads with SPAdes using lower k-mer values and runs in single cell mode.')
        parser.add_argument('--genemark_path', action='store', default=None,
                            help='Path to GeneMark-ES script.', type=str, required=False)
        parser.add_argument('--print_workflow', action='store_true', default="False", required=False, 
                            help='Will print a summary of the workflow upon completion of the script')                                                        

        args = parser.parse_args()

        # Enforce additional argument checks
        if args.singularity_funannotate and not args.funannotate:
            parser.error("argument -sf/--singularity_funannotate: requires -f/--funannotate")
        if args.singularity_antismash and not args.antismash:
            parser.error("argument -sa/--singularity_antismash: requires -ant/--antismash")
        if args.eggnog_db and not args.eggnog:
            parser.error("argument -edb/--eggnog_db: requires -e/--eggnog")
        if args.kraken2_db and not args.kraken2:
            parser.error("argument -kdb/--kraken2_db: requires -k/--kraken")
        if args.its_db and not args.its:
            parser.error("argument -idb/--its_db: requires -its/--its")

    except argparse.ArgumentError:
        lib.print_e("An exception occurred with argument parsing. Check your inputs.")
        parser.print_help(sys.stderr)
        sys.exit()

    return args

### Begin Main Script ###
def main():
    input_args = get_args()
    lib.print_t("⁂⁂⁂⁂⁂⁂⁂⁂ Script Begins ⁂⁂⁂⁂⁂⁂⁂⁂  \n")
    start_time = datetime.datetime.now()
    print(input_args)

    # Sets up class object for filenames
    filenames = lib.Files(input_args)

    # Check databases
    lib.check_databases(input_args, filenames)
    
    ### MAKE INPUT DIR ABSOLUTE PATH
    input_args.directory_new = os.path.abspath(os.path.join(input_args.directory, input_args.array))
    lib.make_path(input_args.directory_new)
    os.chdir(input_args.directory_new)

    # Running 'ASSEMBLY' module
    assembly.main(input_args, filenames)

    # Running 'FUNANNOTATE' module
    if input_args.funannotate is True:
        funannotate.main(input_args, filenames)
    else:
        lib.print_h("Skipping assembly annotation...")

    # Running 'POST_ANALYSIS' module
    post_analysis.main(input_args, filenames)

    lib.print_h(f"Results saved to {filenames.results_csv}")
    lib.print_h("Output files and variables are listed below:")
    filenames.printer()
    print("\n")
    lib.print_tu(f"⁂⁂⁂⁂⁂⁂⁂⁂ Script Finished in {datetime.datetime.now() - start_time} ⁂⁂⁂⁂⁂⁂⁂⁂")


    # TEST CODE FOR GENERATING WORKFLOW DESCRIPTION
    if input_args.print_workflow is True:
        if input_args.kraken2 is True: kraken = " and non-fungal reads removed using kraken2"
        else: kraken = ""
        trimming = f"Short sequence reads were trimmed with Trimmomatic{kraken}."
        if input.args.nanopore is None: 
            trimming + f"MinION reads were trimmed with porechop and length-filtered to {input_args.minimum_length} bp with bbduk.sh before being converted to FASTA format. Short reads were mapped to these reads and corrected using ropebwt2 and FMLRC."
            assembly_proc = "Corrected long reads were assembled with Flye"
        else:
            if input_args.careful is True: parameters = " with the \'careful\' parameter" 
            else: parameters = ""
            assembly_proc = f"Reads were assembled with SPADes{parameters}."
        annotation = f"Annotation was performed using the Funannotate pipeline, using {input_args.min_training_models} for GeneMark."
        if input_args.eggnog is True:
            if input_args.interproscan is True: iprscan = " and InterProScan"
            annotation = annotation + f"Functional annotation was performed using Eggnog{iprscan}."
        post = "Assemblies were evaluated with Quast."
        if input_args.its is True: post = post + "ITS regions were extracted using ITSx."
        if input_args.antismash is True: post = post + "BGCs were predicted using antiSMASH."
        print(f"{trimming} {assembly_proc} {annotation} {post}")

if __name__ == '__main__':
    main()
