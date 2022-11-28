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

    a = f"minimap2 -ax sr {filenames.assembly_fasta} {filenames.trimmedf} {filenames.trimmedr}"
    b = f"samtools sort -@ {input_args.cpus} -o {filenames.bamfile}"
    if len(filenames.singularity) > 0: a = " ".join(filenames.singularity) + " " + a
    if len(filenames.singularity) > 0: b = " ".join(filenames.singularity) + " " + b
    commands = [a,b]
    try:
        lib.execute_shell(commands,stdout,stderr)
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

    cmd = ["samtools","index","-@",input_args.cpus,filenames.bamfile]
    print(" ".join(cmd))
    if len(filenames.singularity) > 0: cmd = filenames.singularity + cmd
    try:
        lib.execute(cmd,stdout,stderr)
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
        
    cmd = f"blastn -task megablast -query {filenames.assembly_fasta} -db {input_args.blob_db} -outfmt \"6 qseqid staxids bitscore std sscinames sskingdoms stitle\" -max_target_seqs 1 -max_hsps 1 -num_threads {input_args.cpus} -evalue 1e-25 -out {filenames.megablast_out}"
    if len(filenames.singularity) > 0: cmd = str(" ".join(filenames.singularity)) + " " + cmd
    commands = [cmd]
    try:
        lib.execute_shell(commands,stdout,stderr)
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

    cmd1 = ["blobtools","create","-i",filenames.assembly_fasta,"-b",filenames.bamfile,"-t",filenames.megablast_out,"-o",os.path.join(blobpath,"assembly")]
    cmd2 = ["blobtools","view","-i",filenames.blob_json,"-o",os.path.join(blobpath,"")]
    cmd3 = ["blobtools","plot","-i",filenames.blob_json,"-o",os.path.join(blobpath,"")]
    if len(filenames.singularity) > 0: cmd1 = filenames.singularity + cmd1
    if len(filenames.singularity) > 0: cmd2 = filenames.singularity + cmd2
    if len(filenames.singularity) > 0: cmd3 = filenames.singularity + cmd3    
    try:
        if lib.file_exists(filenames.blob_json,f"blobtools already prepared blob.json!",f"Need to run blobtools create") is False:
            lib.execute(cmd1,stdout,stderr)
            lib.execute(cmd2,stdout,stderr)
        if lib.file_exists(filenames.blob_png,f"blobtools already plotted!",f"Need to run blobtools plot") is False:
            lib.execute(cmd3,stdout,stderr)
        lib.file_exists_exit(filenames.blob_png,"blobtools completed successfully!","blobtools failed.")
    except subprocess.CalledProcessError as e:
        print(e.returncode)
        print(e.output)

def main(input_args,filenames):

    lib.print_h("Initializing \'blobplot\' module...")
    blob_start_time = datetime.datetime.now()
    if input_args is None:
        input_args = get_args()
        filenames = lib.Files(input_args.array)
        filenames.trimmedf = input_args.trimmed_f
        filenames.trimmedr = input_args.trimmed_r
    
    os.chdir(input_args.directory_new)
    blobplot_text = "  ____________\n___/  Blobplot  \_______________________________________________________________"
    lib.print_t(blobplot_text)
    lib.print_h(f"Generating Blobplot of {filenames.assembly_fasta}")
    lib.print_n(input_args)
    
    blobpath = os.path.join(input_args.directory_new,"blobplots")
    lib.make_path(blobpath)
    filenames.bamfile = os.path.join(blobpath,"assembly.bam")
    filenames.indexfile = os.path.join(blobpath,"assembly.bam.bai")
    filenames.megablast_out = os.path.join(blobpath,"megablast.out")
    filenames.blob_json = os.path.join(blobpath,"assembly.blobDB.json")
    filenames.blob_png = os.path.join(blobpath,"assembly.blobDB.json.bestsum.phylum.p8.span.100.blobplot.read_cov.bam0.png")
    print(f"Blobpath = {blobpath} \nBAM file = {filenames.bamfile}")
    if lib.file_exists(filenames.bamfile, "Trimmed reads already mapped to contigs, skipping...","Mapping trimmed reads to contigs with minimap2") is False:
        minimap2(input_args,filenames,blobpath)
    if lib.file_exists(filenames.indexfile, "BAM file already sorted, skipping...","Sorting BAM file with samtools") is False:
        index(input_args,filenames,blobpath)
    if lib.file_exists(filenames.megablast_out, "Contigs already assigned with megablast, skipping...","Assigning contigs with megablast") is False:
        megablast(input_args,filenames,blobpath)
    if lib.file_exists(filenames.blob_png, "Blobtools already executed, skipping...","Running blobtools on input files") is False:
        blobtools(input_args,filenames,blobpath)
    lib.print_h(f"Blobplot script completed in {datetime.datetime.now() - blob_start_time}")

if __name__ == '__main__':
    main()
