import os
import subprocess
import datetime
import library as lib

"""
Runs a basic set of modules from the Funannotate pipeline, as a module of the Fungiflow pipeline.

An input assembly will be cleaned, sorted, softmasked, then gene conding regions will be predicted, 
returning an output folder containing Funannotate outputs, inlcuding GBK, GFF3, an intermediate files.
For further information about the Funannotate pipeline, please visit:
https://funannotate.readthedocs.io/en/latest/index.html
"""

def funannotate_clean(input_args,filenames):
    """
    Runs `funannotate clean` on assembly file to remove duplicate/short contigs.

    Input:      assembly FASTA file
    Output:     cleaned assembly FASTA file
    """

    stdout = f"{input_args.array}_clean.out"
    stderr = f"{input_args.array}_clean.err"

    lib.print_n("Cleaning assembly with funannotate clean")
    cmd1 = ["singularity","exec","-B","/nfs:/nfs",input_args.singularity_funannotate, \
        "funannotate","clean","-i",filenames.assembly_fasta,"-o",filenames.funannotate_clean_fasta,"--cpus",input_args.cpus]
    try:
        lib.execute(cmd1,stdout,stderr)
        lib.file_exists_exit(filenames.funannotate_clean_fasta, \
            "Assembly successfuly cleaned by Funannotate","Funannotate clean failed... check the logs and your inputs")
    except subprocess.CalledProcessError as e:
        print(e.returncode)
        print(e.output) 

def funannotate_sort(input_args,filenames):
    """
    Runs `funannotate sort` on assembly file to sort contigs by size.

    Input:      cleaned assembly FASTA file
    Output:     sorted and cleaned assembly FASTA file
    """

    stdout = f"{input_args.array}_sort.out"
    stderr = f"{input_args.array}_sort.err"

    lib.print_n("Sorting contigs with funannotate sort")
    cmd1 = ["singularity","exec","-B","/nfs:/nfs",input_args.singularity_funannotate, \
        "funannotate","sort","-i",filenames.funannotate_clean_fasta,"-o",filenames.funannotate_sort_fasta]
    try:
        lib.execute(cmd1,stdout,stderr)
        lib.file_exists_exit(filenames.funannotate_sort_fasta, \
            "Assembly successfuly sorted by Funannotate","Funannotate sort failed... check the logs and your inputs")
    except subprocess.CalledProcessError as e:
        print(e.returncode)
        print(e.output) 

def funannotate_mask(input_args,filenames):
    """
    Runs `funannotate mask` on assembly file to softmask non-genic regions of assembly.

    Input:      sorted and cleaned assembly FASTA file
    Output:     sorted, cleaned and masked assembly FASTA file
    """

    stdout = f"{input_args.array}_mask.out"
    stderr = f"{input_args.array}_mask.err"

    lib.print_n("Soft masking repetitive regions with funannotate mask")
    cmd1 = ["singularity","exec","-B","/nfs:/nfs",input_args.singularity_funannotate, \
        "funannotate","mask","-i",filenames.funannotate_sort_fasta,"-o",filenames.funannotate_mask_fasta,"--cpus",input_args.cpus]
    try:
        lib.execute(cmd1,stdout,stderr)
        lib.file_exists_exit(filenames.funannotate_mask_fasta, \
            "Assembly successfully masked by Funannotate","Funannotate mask failed... check the logs and your inputs")
    except subprocess.CalledProcessError as e:
        print(e.returncode)
        print(e.output) 

def funannotate_predict(input_args,filenames,funannotate_path):
    """
    Runs `funannotate predict` on assembly file to predict genes.

    Input:      sorted, cleaned and masked assembly FASTA file
    Output:     annotated assembly GBK file
    """

    stdout = f"{input_args.array}_predict.out"
    stderr = f"{input_args.array}_predict.err"
    
    lib.print_n("Predicting genes with funannotate predict")
    cmd1 = ["singularity","exec","-B","/nfs:/nfs",input_args.singularity_funannotate, \
        "funannotate","predict","-i",filenames.funannotate_mask_fasta,"-o",funannotate_path,"-s",input_args.array,"--cpus",input_args.cpus]
    try:
        lib.execute(cmd1,stdout,stderr)
        lib.file_exists_exit(filenames.funannotate_gbk, \
            "Assembly genes successfully predicted with Funannotate","Funannotate predict failed... check the logs and your inputs")
    except subprocess.CalledProcessError as e:
        print(e.returncode)
        print(e.output)

def main(input_args,filenames):

    lib.print_h("Initializing \'funannotate\' module...")
    funannotate_start_time = datetime.datetime.now()
    annotation_text = "  ____________\n___/ Annotation \_______________________________________________________________"
    lib.print_t(annotation_text)
    lib.print_h("Annotating assembly with Funannotate")
    funannotate_path = os.path.join(input_args.directory_new,"funannotate")
    lib.make_path(funannotate_path)
    lib.print_n("Changing into Funannotate directory")
    os.chdir(funannotate_path)
    filenames.funannotate_clean_fasta = os.path.join(funannotate_path,f"{input_args.array}_cleaned.fasta")
    filenames.funannotate_sort_fasta = os.path.join(funannotate_path,f"{input_args.array}_sorted.fasta")
    filenames.funannotate_mask_fasta = os.path.join(funannotate_path,f"{input_args.array}_masked.fasta")
    filenames.funannotate_gbk = os.path.join(funannotate_path,"predict_results",f"{input_args.array}.gbk")
    filenames.funannotate_gff = os.path.join(funannotate_path,"predict_results",f"{input_args.array}.gff3")

    if lib.file_exists(filenames.funannotate_clean_fasta,"Assembly already cleaned by Funannotate! Skipping...","") is False:
        funannotate_clean(input_args,filenames)
    if lib.file_exists(filenames.funannotate_sort_fasta,"Assembly already sorted by Funannotate! Skipping...","") is False:
        funannotate_sort(input_args,filenames)
    if lib.file_exists(filenames.funannotate_mask_fasta,"Assembly already masked by Funannotate! Skipping...","") is False:
        funannotate_mask(input_args,filenames)
    if lib.file_exists(filenames.funannotate_gbk,"Assembly genes already predicted with Funannotate! Skipping...","") is False:
        funannotate_predict(input_args,filenames,funannotate_path)
    
    lib.print_n("Changing back to main directory")
    os.chdir(input_args.directory_new)
    lib.print_h(f"funannotate module completed in {datetime.datetime.now() - funannotate_start_time}")

if __name__ == '__main__':
    main()