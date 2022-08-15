#!/usr/bin/env
# -*- coding: utf-8 -*-

import os
import sys, argparse
import datetime
import subprocess
import library as lib

def get_args():
    """Parse command line arguments"""
    
    try:
        parser = argparse.ArgumentParser(
            description="Blobplots - from the Fungiflow pipeline.")
        parser.add_argument('-d', '--directory', action='store',
                            help='Working directory path.', type=str, required=True)     
        parser.add_argument('-if', '--trimmed_f', action='store',
                            help='Trimmed Illumina short forward reads path', type=str, required=True)
        parser.add_argument('-ir', '--trimmed_r', action='store',
                            help='Trimmed Illumina short reverse reads path', type=str, required=True)
        parser.add_argument('-as', '--assembly', action='store',
                            help='Assembly FASTA file', type=str, required=True)
        parser.add_argument('-a', '--array', action='store',
                            help='array value', type=str, required=True)
        parser.add_argument('-c', '--cpus', '--cpus', action='store',
                            help='Number of threads to use', type=str, required=True)
        parser.add_argument('-db', '--database', action='store',
                            help='Path to BLAST nt database', type=str, required=True)
        parser.add_argument('-s', '--singularity', '--singularity', action='store',
                            help='Primary Singularity image for Fungiflow', type=str, required=True)

        if len(sys.argv) < 8:
            parser.print_help(sys.stderr)
            print("Parser expecting 8 arguments")
            exit(1)
    except argparse.ArgumentError:
        print("An exception occurred with argument parsing. Check your inputs.")
        exit(1)

    return parser.parse_args()

def minimap2(input_args,filenames,blobpath):
    """
    Maps short Illumina reads to the assembly file using `minimap2`.

    Input:      assembly FASTA file, forward and reverse trimmed FASTQ read files
    Output:     BAM mapping file
    """

    stdout = os.path.join(blobpath,f"{input_args.array}_map.out")
    stderr = os.path.join(blobpath,f"{input_args.array}_map.err")

    cmd1 = ["singularity","exec",input_args.singularity_image, \
        "minimap2","-ax","sr",filenames.assembly_fasta,filenames.trimmedf,filenames.trimmedr,"|","samtools","sort","-@",input_args.cpus,"-o",filenames.bamfile]
    print(" ".join(cmd1))
    try:
        lib.execute_shell(cmd1,stdout,stderr)
        lib.file_exists_exit(filenames.bamfile,"Minimap2 successfully generated BAM file.","Minimap2 mapping failed.")
    except subprocess.CalledProcessError as e:
        print(e.returncode)
        print(e.output)   

def index(input_args,filenames,blobpath):
    """
    Indexes BAM file using `samtools`. 

    Input:      trimmed reads mapped to assembly file - BAM file
    Output:     BAI index file
    """

    stdout = os.path.join(blobpath,f"{input_args.array}_index.out")
    stderr = os.path.join(blobpath,f"{input_args.array}_index.err")

    cmd1 = ["singularity","exec",input_args.singularity_image, \
        "samtools","index","-@",input_args.cpus,filenames.bamfile,filenames.outfile]
    print(" ".join(cmd1))
    try:
        lib.execute(cmd1,stdout,stderr)
        lib.file_exists_exit(filenames.indexfile,"Samtools successfully indexed BAM file.","Samtools indexing failed.")
    except subprocess.CalledProcessError as e:
        print(e.returncode)
        print(e.output)   

def megablast(input_args,filenames,blobpath):
    """
    Uses `megablast` to search nucleotide database for each contig in an assembly, using an E-value cutoff of 1E-25.
    
    Input:      assembly FASTA file
    Output:     megablast results OUT file
    """

    stdout = os.path.join(blobpath,f"{input_args.array}_blast.out")
    stderr = os.path.join(blobpath,f"{input_args.array}_blast.err")
        
    cmd1 = ["singularity","exec",input_args.singularity_image, \
        "blastn","-task","megablast","-query",filenames.assembly_fasta,"-db",input_args.nt_db,"-outfmt","\"6 qseqid staxids bitscore std\"","-max_target_seqs","1","-max_hsps","1","-num_threads",input_args.cpus,"-evalue","1e-25","-out",filenames.megablast_out]
    try:
        lib.execute(cmd1,stdout,stderr)
        lib.file_exists_exit(filenames.megablast_out,"MegaBLAST successfully identified contigs.","MegaBLAST failed.")
    except subprocess.CalledProcessError as e:
        print(e.returncode)
        print(e.output)

def blobtools(input_args,filenames,blobpath):
    """
    Runs the `blobtools` suite of tools on input files to create blobplots.

    Input:      BAM map file, assembly FASTA file, megablast OUT file
    Output:     blobplot JSON file, blobplot PNG images
    """

    stdout = os.path.join(blobpath,f"{input_args.array}_blob.out")
    stderr = os.path.join(blobpath,f"{input_args.array}_blob.err")

    cmd1 = ["singularity","exec",input_args.singularity_image, \
        "blobtools-blobtools_v1.1.1/blobtools","create","-i",filenames.assembly_fasta,"-b",filenames.bamfile,"-t",filenames.megablast_out,"-o",filenames.blob_json]
    cmd2 = ["singularity","exec",input_args.singularity_image, \
        "blobtools-blobtools_v1.1.1/blobtools","view","-i",filenames.blob_json]
    cmd3 = ["singularity","exec",input_args.singularity_image, \
        "blobtools-blobtools_v1.1.1/blobtools","plot","-i",filenames.blob_json]
    try:
        lib.execute(cmd1,stdout,stderr)
        lib.execute(cmd2,stdout,stderr)
        lib.execute(cmd3,stdout,stderr)
        lib.file_exists_exit(filenames.blob_json,"Blobtools successfully generated files","Blobtools failed.")
    except subprocess.CalledProcessError as e:
        print(e.returncode)
        print(e.output)

def main(input_args,filenames):

    start_time = datetime.datetime.now()
    if input_args is None:
        input_args = get_args()
        filenames = lib.Files(input_args.trimmed_f.split("_")[0])
        filenames.trimmedf = input_args.trimmed_f
        filenames.trimmedr = input_args.trimmed_r
    
    os.chdir(input_args.directory_new)
    blobplot_text = "   ____________\n___/  Blobplot  \_______________________________________________________________"
    lib.print_t(blobplot_text)
    lib.print_h(f"Generating Blobplot of {input_args.assembly}")
    lib.print_n(input_args)
    
    blobpath = os.path.join(input_args.directory_new,"blobplots")
    lib.make_path(blobpath)
    name = input_args.assembly.split(".")[0]
    filenames.bamfile = os.path.join(blobpath,f"{name}.bam")
    filenames.indexfile = os.path.join(blobpath,f"{name}.bai")
    filenames.megablast_out = os.path.join(blobpath,f"{name}.out")
    filenames.blob_json = os.path.join(blobpath,f"{name}.blobDB.json")
    print(f"Blobpath = {blobpath} \nBAM file = {filenames.bamfile}")
    print("Mapping trimmed reads to assembly with Minimap2")
    minimap2(input_args,filenames,blobpath)
    print("Sorting BAM with Samtools")
    index(input_args,filenames,blobpath)
    print("Megablasting assembly")
    megablast(input_args,filenames,blobpath)
    print("Blobbing and plotting")
    blobtools(input_args,filenames,blobpath)
    print(f"Blobplot script completed in {datetime.datetime.now() - start_time}")

if __name__ == '__main__':
    main()