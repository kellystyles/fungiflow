import os
import subprocess
import datetime
import library as lib

"""
Runs a basic set of modules from the Funannotate pipeline, as a module of the Fungiflow pipeline.

    - regarding GeneMark-ES usage, it is not compulsory but does give good predictions.
    If you wish to use this software, please obtain a license key and install it locally,
    and supply the path to `gmes_petap.pl` using the Fungiflow command 

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
    cmd = ["funannotate","clean","-i",filenames.assembly_fasta,"-o",filenames.funannotate_clean_fasta,"--cpus",input_args.cpus]
    if len(filenames.funannotate) > 0: cmd = filenames.funannotate + cmd   
    try:
        lib.execute(cmd,stdout,stderr)
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
    cmd = ["funannotate","sort","-i",filenames.funannotate_clean_fasta,"-o",filenames.funannotate_sort_fasta]
    if len(filenames.funannotate) > 0: cmd = filenames.funannotate + cmd
    try:
        lib.execute(cmd,stdout,stderr)
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
    cmd = ["funannotate","mask","-i",filenames.funannotate_sort_fasta,"-o",filenames.funannotate_mask_fasta,"--cpus",input_args.cpus]
    if len(filenames.funannotate) > 0: cmd = filenames.funannotate + cmd
    try:
        lib.execute(cmd,stdout,stderr)
        lib.file_exists_exit(filenames.funannotate_mask_fasta, \
            "Assembly successfully masked by Funannotate","Funannotate mask failed... check the logs and your inputs")
    except subprocess.CalledProcessError as e:
        print(e.returncode)
        print(e.output) 

def funannotate_predict(input_args,filenames,funannotate_path):
    """
    Runs `funannotate predict` on assembly file to predict genes.
    If GeneMark-ES doesn't have enough training models, causing Funannotate to fail,
    Will run GeneMark-ES standalone with parameters optimized for smaller contig
    sizes and feed the output file to funannotate predict again.

    Input:      sorted, cleaned and masked assembly FASTA file
    Output:     annotated assembly GBK file
    """

    stdout = f"{input_args.array}_predict.out"
    stderr = f"{input_args.array}_predict.err"
    
    lib.print_n("Predicting genes with funannotate predict")
    cmd = ["funannotate","predict","-i",filenames.funannotate_mask_fasta, \
        "-o",funannotate_path,"-s",input_args.array,"--cpus",input_args.cpus, \
            "--name",input_args.array,"--min_training_models",input_args.min_training_models]
    if len(filenames.funannotate) > 0: cmd = filenames.funannotate + cmd
    if input_args.genemark_path is not None: cmd = cmd + ["--GENEMARK_PATH", input_args.genemark_path]
    try:
        lib.execute(cmd,stdout,stderr)
        # If funannotate fails because GeneMark-ES doesn't have enough training models,
        # run GeneMark-ES separately with reduced contig sizes
        if lib.file_exists(filenames.funannotate_gbk, \
            "Assembly genes successfully predicted with Funannotate predict","") is False:
            with open(stderr, "r") as f:
                lines = f.readlines()
                if "error, file not found: data/training.fna" in lines[-2]:
                    gme_cmd = ["gmes_petap.pl","--ES","--max_intron","3000","--soft_mask","2000","--cores",input_args.cpus,"--sequence", \
                        os.path.join(funannotate_path,"predict_misc","genome.softmasked.fa"),"--fungus","--min_contig","2000", \
                            "--work_dir",os.path.join(funannotate_path,"predict_misc")]
                    if len(filenames.funannotate) > 0: gme_cmd = filenames.funannotate + gme_cmd
                    lib.print_n("GeneMark ES could not predict genes from the assembly as it is too fragmented. Will try to optimize the parameters for low-coverage assembly...")
                    lib.execute(gme_cmd,stdout,stderr)
                    cmd = cmd + ["--genemark_gtf",os.path.join(funannotate_path,"predict_misc","genemark.gtf")]
                    lib.execute(cmd,stdout,stderr)
                lib.file_exists(filenames.funannotate_gbk, \
                    "Assembly genes successfully predicted with Funannotate predict", \
                        "Funannotate predict failed... check the logs and your inputs\nYou could try using the `--careful` input parameter to prepare a better assembly. You may need to remove the existing `assembly` directory first.")

    except subprocess.CalledProcessError as e:
        print(e.returncode)
        print(e.output)

def eggnog_annotate(input_args,filenames,eggnog_path):
    """
    Runs eggnog mapper to functionally annotate predicted protein sequences.

    Input:      predicted proteins multi FASTA
    Output:     eggnog annotations of proteins multi FASTA
    """

    stdout = os.path.join(eggnog_path,f"{input_args.array}_eggnog_annotate.out")
    stderr = os.path.join(eggnog_path,f"{input_args.array}_eggnog_annotate.err")
    output_prefix = os.path.join(eggnog_path,input_args.array)
    lib.print_n("Annotating CDS with eggnog mapper")
    cmd = ["emapper.py","-m","diamond","-i",filenames.funannotate_prots,"--data_dir",input_args.eggnog_db, \
        "-o",output_prefix,"--tax_scope","4751,33154,2759,1"]
    if len(filenames.funannotate) > 0: cmd = filenames.funannotate + cmd
    # If there is sufficient memory, load the entire diamond DB into memory
    if int(input_args.mem) > 45: cmd = cmd + ["--dbmem"]
    # If a hits file from a previous run is present, resume the run
    if os.path.isfile(os.path.join(eggnog_path, f"{input_args.array}.emapper.hits")) is True: cmd = cmd + ["--resume"]
    try:
        lib.execute(cmd,stdout,stderr)
        lib.file_exists(filenames.eggnog_annotations, \
            "CDS successfully annotated with eggnog mapper","eggnog mapper failed... check the logs and your inputs")
    except subprocess.CalledProcessError as e:
        print(e.returncode)
        print(e.output)

def interproscan(input_args,filenames,iprscan_path):
    """
    Runs InterProScan to functionally annotate the GFF file.

    Input: GFF3 file - output of funannotate predict
    Output: XML file with functional hits
    """

    stdout = f"{input_args.array}_interproscan.out"
    stderr = f"{input_args.array}_interproscan.err"
    lib.print_n("Annotating CDS with InterProScan")
    cmd = ["interproscan.sh","-i",input_args.funannotate_prots, \
        "-d",iprscan_path,"-f","XML","-goterms","-iprlookup","-pa","-cpu",input_args.cpus]
    try:
        lib.execute(cmd,stdout,stderr)
        lib.file_exists(filenames.funannotate_func_gbk, \
            "CDS successfully annotated with Funannotate annotate", \
                "Funannotate annotate failed... check the logs and your inputs")
    except subprocess.CalledProcessError as e:
        print(e.returncode)
        print(e.output)
   
def funannotate_annotate(input_args,filenames,funannotate_path):
    """
    Runs `funannotate predict` on assembly file to predict genes.

    Input:      Funannotate predict_results folder
    Output:     functionally annotated assembly GBK file
    """

    stdout = f"{input_args.array}_annotate.out"
    stderr = f"{input_args.array}_annotate.err"
    
    lib.print_n("Annotating CDS with funannotate annotate")
    cmd = ["funannotate","annotate","-i",os.path.join(funannotate_path,"predict_results"), \
        "-o",funannotate_path,"-s",input_args.array,"--cpus",input_args.cpus]
    if lib.file_exists(filenames.eggnog_annotations,"","") is True: 
        cmd = cmd + ["--eggnog",filenames.eggnog_annotations]
        try:
            os.remove(os.path.join(funannotate_path, "annotate_misc", "eggnog.emapper.annotations"))
        except FileNotFoundError:
            pass
    if len(filenames.funannotate) > 0: cmd = filenames.funannotate + cmd
    try:
        lib.execute(cmd,stdout,stderr)
        lib.file_exists(filenames.funannotate_func_gbk, \
            "CDS successfully annotated with Funannotate annotate", \
                "Funannotate annotate failed... check the logs and your inputs")
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
    eggnog_path = os.path.join(funannotate_path, "eggnog")
    iprscan_path = os.path.join(funannotate_path, "interproscan")
    lib.make_path(funannotate_path)
    lib.make_path(eggnog_path)
    lib.print_n("Changing into Funannotate directory")
    os.chdir(funannotate_path)

    # Defining filenames
    filenames.funannotate_clean_fasta = os.path.join(funannotate_path,f"{input_args.array}_cleaned.fasta")
    filenames.funannotate_sort_fasta = os.path.join(funannotate_path,f"{input_args.array}_sorted.fasta")
    filenames.funannotate_mask_fasta = os.path.join(funannotate_path,f"{input_args.array}_masked.fasta")
    filenames.funannotate_gbk = os.path.join(funannotate_path,"predict_results",f"{input_args.array}.gbk")
    filenames.funannotate_gff = os.path.join(funannotate_path,"predict_results",f"{input_args.array}.gff3")
    filenames.funannotate_prots = os.path.join(funannotate_path,"predict_results",f"{input_args.array}.proteins.fa")
    filenames.eggnog_annotations = os.path.join(eggnog_path,f"{input_args.array}.emapper.annotations")
    filenames.ipr_annotations = os.path.join(iprscan_path,f"{input_args.array}.xml")
    filenames.funannotate_func_gbk = os.path.join(funannotate_path,"annotate_results",f"{input_args.array}.gbk")

    # Check if output file exists or not before running each command
    if lib.file_exists(filenames.funannotate_clean_fasta,"Assembly already cleaned by Funannotate! Skipping...","") is False:
        funannotate_clean(input_args,filenames)
    if lib.file_exists(filenames.funannotate_sort_fasta,"Assembly already sorted by Funannotate! Skipping...","") is False:
        funannotate_sort(input_args,filenames)
    if lib.file_exists(filenames.funannotate_mask_fasta,"Assembly already masked by Funannotate! Skipping...","") is False:
        funannotate_mask(input_args,filenames)
    if lib.file_exists(filenames.funannotate_gbk,"Assembly genes already predicted with Funannotate! Skipping...","") is False:
        funannotate_predict(input_args,filenames,funannotate_path)
    if lib.file_exists(filenames.funannotate_gbk,"","") is True \
        and lib.file_exists(filenames.eggnog_annotations,"Functional annotation already completed with eggnog","") is False:
        eggnog_annotate(input_args,filenames,eggnog_path)
    if lib.file_exists(filenames.funannotate_func_gbk,"Functional annotation with Funannotate! Skipping...","") is False:
        funannotate_annotate(input_args,filenames,funannotate_path)

    lib.print_n("Changing back to main directory")
    os.chdir(input_args.directory_new)
    lib.print_h(f"funannotate module completed in {datetime.datetime.now() - funannotate_start_time}")

if __name__ == '__main__':
    main()