import os
import subprocess
import library as lib

"""
Runs a basic set of modules from the Funannotate pipeline, as a module of the Fungiflow pipeline.

An input assembly will be cleaned, sorted, softmasked, then gene conding regions will be predicted, 
returning an output folder containing Funannotate outputs, inlcuding GBK, GFF3, an intermediate files.
For further information about the Funannotate pipeline, please visit:
https://funannotate.readthedocs.io/en/latest/index.html
"""

def funannotate_clean(input_args,assembly,funannotate_clean):
    """
    Runs `funannotate clean` on assembly file to remove duplicate/short contigs.

    Input:      assembly FASTA file
    Output:     cleaned assembly FASTA file
    """

    stdout = f"{input_args.array}_clean.out"
    stderr = f"{input_args.array}_clean.err"

    lib.print_n(" - cleaning assembly")
    cmd1 = ["singularity","exec","-B","/nfs:/nfs",input_args.singularity_funannotate, \
        "funannotate","clean","-i",assembly,"-o",funannotate_clean,"--cpus",input_args.cpus]
    try:
        lib.execute(cmd1,stdout,stderr)
        lib.file_exists_exit(funannotate_clean_fasta, \
            "Assembly successfuly cleaned by Funannotate","Funannotate clean failed... check the logs and your inputs")
    except subprocess.CalledProcessError as e:
        print(e.returncode)
        print(e.output) 

def funannotate_sort(input_args,funannotate_clean,funannotate_sort):
    """
    Runs `funannotate sort` on assembly file to sort contigs by size.

    Input:      cleaned assembly FASTA file
    Output:     sorted and cleaned assembly FASTA file
    """

    stdout = f"{input_args.array}_sort.out"
    stderr = f"{input_args.array}_sort.err"

    lib.print_n(" - sorting contigs")
    cmd1 = ["singularity","exec","-B","/nfs:/nfs",input_args.singularity_funannotate, \
        "funannotate","sort","-i",funannotate_clean,"-o",funannotate_sort]
    try:
        lib.execute(cmd1,stdout,stderr)
        lib.file_exists_exit(funannotate_sort_fasta, \
            "Assembly successfuly sorted by Funannotate","Funannotate sort failed... check the logs and your inputs")
    except subprocess.CalledProcessError as e:
        print(e.returncode)
        print(e.output) 

def funannotate_mask(input_args,funannotate_sort,funannotate_mask):
    """
    Runs `funannotate mask` on assembly file to softmask non-genic regions of assembly.

    Input:      sorted and cleaned assembly FASTA file
    Output:     sorted, cleaned and masked assembly FASTA file
    """

    stdout = f"{input_args.array}_mask.out"
    stderr = f"{input_args.array}_mask.err"

    lib.print_n(" - soft masking repetitive regions")
    cmd1 = ["singularity","exec","-B","/nfs:/nfs",input_args.singularity_funannotate, \
        "funannotate","mask","-i",funannotate_sort,"-o",funannotate_mask,"--cpus",input_args.cpus]
    try:
        lib.execute(cmd1,stdout,stderr)
        lib.file_exists_exit(funannotate_mask_fasta, \
            "Assembly successfuly masked by Funannotate","Funannotate mask failed... check the logs and your inputs")
    except subprocess.CalledProcessError as e:
        print(e.returncode)
        print(e.output) 

def funannotate_predict(input_args,funannotate_path,funannotate_mask):
    """
    Runs `funannotate predict` on assembly file to predict genes.

    Input:      sorted, cleaned and masked assembly FASTA file
    Output:     annotated assembly GBK file
    """

    stdout = f"{input_args.array}_predict.out"
    stderr = f"{input_args.array}_predict.err"

    cmd1 = ["singularity","exec","-B","/nfs:/nfs",input_args.singularity_funannotate, \
        "funannotate","predict","-i",funannotate_mask,"-o",funannotate_path,"-s",input_args.array,"--cpus",input_args.cpus]
    try:
        lib.execute(cmd1,stdout,stderr)
        lib.file_exists_exit(funannotate_gbk, \
            "Assembly genes successfully predicted with Funannotate","Funannotate predict failed... check the logs and your inputs")
    except subprocess.CalledProcessError as e:
        print(e.returncode)
        print(e.output)

def main(input_args,filenames):

    annotation_text = "   ____________\n___/ Annotation \_______________________________________________________________"
    lib.print_t(annotation_text)
    lib.print_h("Annotating assembly with Funannotate")
    funannotate_path = os.path.join(directory_new,"funannotate")
    lib.make_path(funannotate_path)
    lib.print_n("Changing into Funannotate directory")
    os.chdir(funannotate_path)
    filenames.funannotate_clean_fasta = os.path.join(funannotate_path,f"{input_args.array}_cleaned.fasta")
    filenames.funannotate_sort_fasta = os.path.join(funannotate_path,f"{input_args.array}_sorted.fasta")
    filenames.funannotate_mask_fasta = os.path.join(funannotate_path,f"{input_args.array}_masked.fasta")
    filenames.funannotate_gbk = os.path.join(funannotate_path,"predict_results",f"{input_args.array}.gbk")
    filenames.funannotate_gff = os.path.join(funannotate_path,"predict_results",f"{input_args.array}.gff3")

    if lib.file_exists(funannotate_clean_fasta,"Assembly already cleaned by Funannotate! Skipping...","Cleaning assembly with Funannotate clean") is False:
        funannotate_clean(input_args,filenames.assembly_fasta,funannotate_clean_fasta)
    if lib.file_exists(funannotate_sort_fasta,"Assembly already sorted by Funannotate! Skipping...","Sorting assembly with Funannotate sort") is False:
        funannotate_sort(input_args,funannotate_clean_fasta,funannotate_sort_fasta)
    if lib.file_exists(funannotate_mask_fasta,"Assembly already masked by Funannotate! Skipping...","Masking repetitive regions with Funannotate mask") is False:
        funannotate_mask(input_args,funannotate_sort_fasta,funannotate_mask_fasta)
    if lib.file_exists(funannotate_gbk,"Assembly genes already predicted with Funannotate! Skipping...","Predicting genes with Funannotate predict") is False:
        funannotate_predict(input_args,funannotate_path,funannotate_mask_fasta)
    
    lib.print_n("Changing back to main directory")
    os.chdir(input_args.directory_new)
  
    