import os
import subprocess
import datetime
import library as lib

"""
The assembly module of the Fungiflow pipeline.
"""

def find_adapters(input_args,filenames,adapter_path):
    """
    Finds adapter sequences in raw read files.  
    
    Inputs:     forward and reverse raw Illumina paired-end FASTQ reads.
    Outputs:    FASTA file containing identified adapter sequences.
    """

    stdout = os.path.join(adapter_path,f"adapter_{input_args.array}.out")
    stderr = os.path.join(adapter_path,f"adapter_{input_args.array}.err")

    lib.print_h(f"Identifying adapter sequences in {input_args.array} reads")
    cmd1 = ["singularity","exec","-B","/nfs:/nfs",input_args.singularity,"bash","adap_ID.sh",filenames.shortf,filenames.adapter_file]
    cmd2 = ["singularity","exec","-B","/nfs:/nfs",input_args.singularity,"bash","adap_ID.sh",filenames.shortr,filenames.adapter_file]
    try:
        print(" ".join(cmd1))
        lib.execute(cmd1,stdout,stderr)
        print(" ".join(cmd2))
        lib.execute(cmd2,stdout,stderr)

    except subprocess.CalledProcessError as e:
        print(e.returncode)
        print(e.output)   

def trimmomatic(input_args,filenames,trimmed_path):
    """
    Trims primer sequences from raw Illumina reads with `trimmomatic`.  
    
    Inputs:     forward and reverse Illumina paired-end FASTQ reads.
    Outputs:    forward and reverse trimmed  Illumina paired-end FASTQ reads.
    """

    stdout = os.path.join(trimmed_path,f"{input_args.array}.out")
    stderr = os.path.join(trimmed_path,f"{input_args.array}.err")

    lib.print_h(f"Trimming Illumina {input_args.array} reads with Trimmomatic")
    cmd1 = ["singularity","exec","-B","/nfs:/nfs",input_args.singularity,"TrimmomaticPE","-threads",input_args.cpus,"-trimlog",filenames.trimlog_file,filenames.shortf,filenames.shortr,filenames.trimmedf,filenames.utrimmedf,filenames.trimmedr,filenames.utrimmedr,f"ILLUMINACLIP:{filenames.adapter_file}:2:30:10:4:4:/true","TRAILING:9","SLIDINGWINDOW:4:15","MINLEN:36"]
    print(" ".join(cmd1))
    try:
        lib.execute(cmd1,stdout,stderr)
    except subprocess.CalledProcessError as e:
        print(e.returncode)
        print(e.output)

def fastqc(input_args,trimmed_path):
    """
    Performs QC analysis of Illumina reads with `fastqc`.  
    
    Inputs:     forward and reverse Illumina paired-end FASTQ reads.
    Outputs:    fastqc output folder
    """

    stdout = os.path.join(trimmed_path,f"{input_args.array}_fastqc.out")
    stderr = os.path.join(trimmed_path,f"{input_args.array}_fastqc.err")

    lib.print_n("Assessing reads with FastQC")
    cmd1 = ["singularity","exec","-B","/nfs:/nfs",input_args.singularity,"fastqc","-o",trimmed_path,os.path.join(trimmed_path,f"{input_args.array}_*.fq*")]
    print(" ".join(cmd1))
    try:
        lib.execute(cmd1,stdout,stderr)
    except subprocess.CalledProcessError as e:
        print(e.returncode)
        print(e.output)   

def porechop(input_args,filenames,trimmed_path):
    stdout = os.path.join(trimmed_path,f"{input_args.array}_porechop.out")
    stderr = os.path.join(trimmed_path,f"{input_args.array}_porechop.err")
    # Trimming long reads with PORECHOP
    cmd1 = ["singularity","exec","-B","/nfs:/nfs",input_args.singularity,"porechop","-i",filenames.nanopore,"-o",filenames.nanopore_trimmed,"--adapter_threshold","96","--no_split"]
    try:
        lib.execute(cmd1,stdout,stderr)
    except subprocess.CalledProcessError as e:
        print(e.returncode)
        print(e.output) 

def fastq_2_fasta(input_args,filenames,trimmed_path):
    stdout = os.path.join(trimmed_path,f"{input_args.array}_fastq_2_fasta.out")
    stderr = os.path.join(trimmed_path,f"{input_args.array}_fastq_2_fasta.err")
    # Converting long reads to FASTA format
    cmd1 = ["singularity","exec","-B","/nfs:/nfs",input_args.singularity,"reformat.sh","-f",f"in={filenames.nanopore_trimmed}",f"out={filenames.nanopore_fasta}","ignorebadquality"] 
    try:
        lib.execute(cmd1,stdout,stderr)
    except subprocess.CalledProcessError as e:
        print(e.returncode)
        print(e.output) 

def length_filter_fasta(input_args,filenames,trimmed_path):
    stdout = os.path.join(trimmed_path,f"{input_args.array}_length_filter_fasta.out")
    stderr = os.path.join(trimmed_path,f"{input_args.array}_length_filter_fasta.err")
    # Length-filtering trimmed long reads
    cmd1 = ["singularity","exec","-B","/nfs:/nfs",input_args.singularity,"bbduk.sh","-f",f"in={filenames.nanopore_fasta}",f"out={filenames.nanopore_length_filtered}",f"minlen={input_args.minlen}"]
    try:
        lib.execute(cmd1,stdout,stderr)
    except subprocess.CalledProcessError as e:
        print(e.returncode)
        print(e.output) 

def short_map_2_long(input_args,filenames,trimmed_path):
    stdout = os.path.join(trimmed_path,f"{input_args.array}_short_map_2_long.out")
    stderr = os.path.join(trimmed_path,f"{input_args.array}_short_map_2_long.err")
    # mapping short reads
    cmd1 = ["gunzip","-c",filenames.trimmedf,filenames.trimmedr,"|","awk","\'NR % 4 == 2\'","|","sort","|","tr","NT","TN","|","singularity","exec","-B","/nfs:/nfs",input_args.singularity,"ropebwt2","-LR","|","tr","NT","TN","|","singularity","exec","-B","/nfs:/nfs",input_args.singularity,"fmlrc-convert",filenames.short_mapped_2_long]
    try:
        lib.execute(cmd1,stdout,stderr)
    except subprocess.CalledProcessError as e:
        print(e.returncode)
        print(e.output) 

def long_read_corr(input_args,filenames,trimmed_path):
    stdout = os.path.join(trimmed_path,f"{input_args.array}_long_read_corr.out")
    stderr = os.path.join(trimmed_path,f"{input_args.array}_long_read_corr.err")
    # Correcting long reads with FMLRC
    cmd1 = ["singularity","exec","-B","/nfs:/nfs",input_args.singularity,"fmlrc","-p",input_args.cpus,filenames.short_mapped_2_long,filenames.nanopore_fasta,filenames.nanopore_corr]
    try:
        lib.execute(cmd1,stdout,stderr)
    except subprocess.CalledProcessError as e:
        print(e.returncode)
        print(e.output) 

def long_read_trim(input_args,filenames,trimmed_path):
    """
    Prepares long reads for assembly by:
    1. trim long reads with `porechop`
    2. length-filtering reads, removing reads shorter than 3000 bp with `bbduk.sh`
    3. converting the long reads FASTQ file to FASTA file with `reformat.sh`
    4. mapping short paired Illumina reads to the long reads using `ropebwt2` and `fmlrc-convert`
    5. correcting the long reads with mapped Illumina reads using `fmlrc`
    
    Input:      raw Nanopore FASTQ reads file.
    Output:     corrected Nanopore FASTA reads file.
    """
    stdout = os.path.join(trimmed_path,f"{input_args.array}_long_read_trim.out")
    stderr = os.path.join(trimmed_path,f"{input_args.array}_long_read_trim.err")

    lib.print_h(f"Preparing long reads {input_args.nanopore} for assembly")
    porechop(input_args,filenames,trimmed_path)
    lib.file_exists_exit(filenames.nanopore_trimmed,"Porechop successfully trimmed the reads","Porechop did not complete successfully. Check the logs.")
    fastq_2_fasta(input_args,filenames,trimmed_path)
    lib.file_exists_exit(filenames.nanopore_fasta,"reformat.sh successfully converted the reads to FASTA format","reformat.sh did not complete successfully. Check the logs.")
    length_filter_fasta(input_args,filenames,trimmed_path)
    lib.file_exists_exit(filenames.nanopore_length_filtered,"bbduk.sh successfully length-filtered the reads","bbduk.sh did not complete successfully. Check the logs.")
    short_map_2_long(input_args,filenames,trimmed_path)
    lib.file_exists_exit(filenames.short_mapped_2_long,"The short reads were successfully mapped to the long reads","The short read mapping to long reads was not successful. Check the logs.")
    long_read_corr(input_args,filenames,trimmed_path)
    lib.file_exists_exit(filenames.nanopore_corr,"FMLRC successfully corrected the reads","FMLRC did not complete successfully. Check the logs.")

def assembly_short(input_args,filenames,assembly_path):
    """
    Performs assembly of short reads with `spades` assembler.  
    
    Inputs:     forward and reverse trimmed Illumina paired-end FASTQ reads.
    Outputs:    spades output folder
    """

    stdout = os.path.join(assembly_path,f"{input_args.array}.out")
    stderr = os.path.join(assembly_path,f"{input_args.array}.err")

    if input_args.data_type is "metagenomic":
        assembly_type = "--meta"
    else:
        assembly_type = "--isolate"

    # Assembly with SPADES
    lib.print_h(f"Assembling trimmed {input_args.array} reads with SPADES")
    cmd1 = ["singularity","exec","-B","/nfs:/nfs",input_args.singularity,"spades.py",assembly_type,"--threads",input_args.cpus,"--memory",input_args.mem,"-k","21,33,55,77,99,127","--pe1-1",filenames.trimmedf,"--pe1-2",filenames.trimmedr,"-o",assembly_path]
    print(" ".join(cmd1))
    try:
        lib.execute(cmd1,stdout,stderr)
    except subprocess.CalledProcessError as e:
        print(e.returncode)
        print(e.output)   

def assembly_hybrid(input_args,filenames,assembly_path):
    """
    Performs a hybrid assembly using `flye` assembler.

    Inputs:     forward and reverse trimmed Illumina paired-end FASTQ reads, trimmed Nanopore FASTA reads.
    Outputs:    flye output folder
    """
    stdout = os.path.join(assembly_path,f"{input_args.array}.out")
    stderr = os.path.join(assembly_path,f"{input_args.array}.err")

    if input_args.type is "metagenomic":
        assembly_type = "--meta"
    else:
        assembly_type = ""

    # Assembly with FLYE
    lib.print_h(f"Hybrid isolate assembly of {input_args.array} reads with FLYE")
    cmd1 = ["singularity","exec","-B","/nfs:/nfs",input_args.singularity,"flye","--nano-raw",filenames.nanopore_corr,"-o",assembly_path,"-t",input_args.cpus,"-i","2",assembly_type]
    print(" ".join(cmd1))
    try:
        lib.execute(cmd1,stdout,stderr)
        print(f"{datetime.datetime.now():%Y-%m-%d %I:%M:%S} Completed assembly of trimmed {input_args.array} reads with FLYE")
    except subprocess.CalledProcessError as e:
        print(e.returncode)
        print(e.output)   

def hybrid_polish(input_args,filenames,assembly_path):
    """
    Polishes hybrid assemblies as follows:
    1. map trimmed/length-filtered long reads to assembly with `minimap2`
    2. one round of polishing using `racon`
    3. one round of polishing using `medaka`
    4. map short trimmed/taxa-filtered short reads to medaka output assembly
    5. one round of polishing using `pilon`
    
    Inputs:    hybrid flye assembly FASTA file, Nanopore FASTA and Illumina FASTQ reads.
    Output:    polished assembly FASTA file.
    """
    stdout = os.path.join(assembly_path,f"{input_args.array}.out")
    stderr = os.path.join(assembly_path,f"{input_args.array}.err")
    
    # Polishing FLYE assembly
    lib.print_h(f"Polishing {filenames.assembly_fasta}")
    cmd1 = ["singularity","exec","-B","/nfs:/nfs",input_args.singularity,"minimap2","-x","map-ont",filenames.assembly_fasta,filenames.nanopore_length_filtered,">",filenames.nanopore_sam]
    racon_path = os.path.join(assembly_path,"racon")
    lib.make_path(racon_path)
    os.chdir(racon_path)
    cmd2 = ["singularity","exec","-B","/nfs:/nfs",input_args.singularity,"racon","--no-trimming","-t",input_args.cpus,filenames.nanopore_length_filtered,filenames.nanopore_sam,filenames.assembly_fasta]
    os.chdir(input_args.directory)
    medaka_path = os.path.join(assembly_path,"medaka")
    lib.make_path(medaka_path)
    os.chdir(medaka_path)
    cmd3 = ["singularity","exec","-B","/nfs:/nfs",input_args.singularity,"medaka_consenus","-i",filenames.nanopore_length_filtered,"-d",filenames.racon_assembly,"-o","medaka","-t",input_args.cpus,"-m","r941_min_hac_g511"]
    os.chdir(input_args.directory)
    cmd4 = ["singularity","exec","-B","/nfs:/nfs",input_args.singularity,"bwa-mem2","index",filenames.assembly_fasta]
    cmd4 = ["singularity","exec","-B","/nfs:/nfs",input_args.singularity,"bwa-mem2","mem","-t",input_args.cpus,filenames.medaka_assembly,filenames.trimmedf,filenames.trimmedr,">",filenames.assembly_bam]
    cmd5 = ["singularity","exec","-B","/nfs:/nfs",input_args.singularity,"samtools","sort","-@",input_args.cpus,"-o",filenames.assembly_sorted_bam,filenames.assembly_bam]
    cmd6 = ["singularity","exec","-B","/nfs:/nfs",input_args.singularity,"samtools","index","-@",input_args.cpus,filenames.assembly_sorted_bam]
    cmd7 = ["singularity","exec","-B","/nfs:/nfs",input_args.singularity,"java",f"-Xmx{input_args.mem}G","-jar","pilon","--genome",filenames.assembly_fasta,"--bam",filenames.assembly_sorted_bam,"--outdir",assembly_path]

    print(" ".join(cmd1))
    print(" ".join(cmd2))
    print(" ".join(cmd3))
    print(" ".join(cmd4))
    print(" ".join(cmd5))
    print(" ".join(cmd6))
    print(" ".join(cmd7))
    try:
        lib.execute(cmd1,stdout,stderr)
        lib.execute(cmd2,stdout,stderr)
        lib.execute(cmd3,stdout,stderr)
        lib.execute(cmd4,stdout,stderr)
        lib.execute(cmd5,stdout,stderr)
        lib.execute(cmd6,stdout,stderr)
        lib.execute(cmd7,stdout,stderr)
    except subprocess.CalledProcessError as e:
        print(e.returncode)
        print(e.output)

def kraken2(input_args,filenames,kraken_path):
    """
    Taxonomically filters trimmed Illumina reads using `kraken2`.

    Inputs:     forward and reverse trimmed Illumina reads.
    Outputs:    forward and reverse filtered and unfiltered Illumina reads.
    """
    
    stdout = os.path.join(kraken_path,f"{input_args.array}.out")
    stderr = os.path.join(kraken_path,f"{input_args.array}.err")

    lib.print_h(f"Performing Kraken2 analysis of trimmed reads for {input_args.array}")
    cmd1 = ["singularity","exec","-B","/nfs:/nfs",input_args.singularity,"kraken2","--db",input_args.kraken2_db,\
        "--threads",input_args.cpus,"--use-names","--report",filenames.kreport,"--classified-out",filenames.k_classified,\
            "--unclassified-out",filenames.k_unclassified,"--paired",filenames.trimmedf,filenames.trimmedr]
    print(" ".join(cmd1))    
    try:
        lib.execute(cmd1,stdout,stderr)
    except subprocess.CalledProcessError as e:
        print(e.returncode)
        print(e.output)   

def main(input_args,filenames):
    
    trim_text = "   __________\n___/ Trimming \_________________________________________________________________"
    lib.print_t(trim_text)

    adapter_path = "adapters"
    lib.make_path(adapter_path)

    filenames.adapter_file = os.path.join("adapters",f"{str(input_args.array)}_adapters.fasta")
    if lib.file_exists(filenames.adapter_file,"Adapters already identified! Skipping...","Searching for adapter sequences in short reads") is False:
        find_adapters(input_args,filenames,adapter_path)
        lib.file_exists_exit(filenames.adapter_file,"Adapters identified!","Adapters not successfully extracted... check the logs and your inputs")

    trimmed_path = "trimmed"
    lib.make_path(trimmed_path)

    filenames.trimmedf = os.path.join(trimmed_path,f"{input_args.array}_trimmed_1P.fq.gz")
    filenames.trimmedr = os.path.join(trimmed_path,f"{input_args.array}_trimmed_2P.fq.gz")
    filenames.utrimmedf = os.path.join(trimmed_path,f"{input_args.array}_trimmed_1U.fq.gz")
    filenames.utrimmedr = os.path.join(trimmed_path,f"{input_args.array}_trimmed_2U.fq.gz")
    filenames.trimlog_file = os.path.join(trimmed_path,f"{input_args.array}_trimmomatic.log")

    if lib.file_exists_bool(filenames.trimmedf,filenames.shortf,0.6,"Forward reads already trimmed! Skipping...","") is False and \
        lib.file_exists_bool(filenames.trimmedr,filenames.shortr,0.6,"Reverse reads already trimmed! Skipping...","") is False:
        trimmomatic(input_args,filenames,trimmed_path)
        lib.file_exists_exit(filenames.trimmedf,"Forward reads trimmed successfully!","Forward reads not trimmed... check the logs and your inputs")
        lib.file_exists_exit(filenames.trimmedr,"Reverse reads trimmed successfully!","Reverse reads not trimmed... check the logs and your inputs")

    if input_args.nanopore is not None:
        filenames.nanopore = input_args.nanopore
        filenames.nanopore_trimmed = os.path.join(trimmed_path,filenames.nanopore.split(".")[0]+".fq")
        filenames.nanopore_fasta = os.path.join(trimmed_path,filenames.nanopore.split(".")[0]+".fa")
        filenames.nanopore_length_filtered = os.path.join(trimmed_path,filenames.nanopore.split(".")[0]+"_lf.fa")
        filenames.short_mapped_2_long = os.path.join(trimmed_path,input_args.array+"_mapped.npy")
        filenames.nanopore_corr = os.path.join(trimmed_path,filenames.nanopore.split(".")[0]+"_corr.fa")
        if lib.file_exists(filenames.nanopore_corr,"Long reads already trimmed and corrected! Skipping...","") is False:
            long_read_trim(input_args,filenames,trimmed_path)
            lib.file_exists_exit(filenames.nanopore_corr,"Long reads trimmed and corrected!","Long reads could not be trimmed and corrected. Please check the logs")

    fastqc(input_args,trimmed_path)

    if len(input_args.kraken2_db) > 0:

        tax_filter_text = "   _____________________\n___/ Taxonomic Filtering \______________________________________________________"
        lib.print_t(tax_filter_text)

        kraken2_path = "kraken2"
        lib.make_path(kraken2_path)
        filenames.kreport = os.path.join(kraken2_path,f"{input_args.array}.report")
        filenames.k_classified = os.path.join(kraken2_path,f"{input_args.array}_class#.fq")
        filenames.k_unclassified = os.path.join(kraken2_path,f"{input_args.array}_unclass#.fq")
        filenames.k_classifiedf = os.path.join(kraken2_path,f"{input_args.array}_class_1.fq")
        filenames.k_classifiedr = os.path.join(kraken2_path,f"{input_args.array}_class_2.fq")
        filenames.k_unclassifiedf = os.path.join(kraken2_path,f"{input_args.array}_unclass_1.fq")
        filenames.k_unclassifiedr = os.path.join(kraken2_path,f"{input_args.array}_unclass_2.fq")

        if lib.file_exists_bool(filenames.k_unclassifiedf,filenames.trimmedf,0.6,"Trimmed reads already filtered by Kraken2! Skipping...","") is False \
        and lib.file_exists_bool(filenames.k_unclassifiedr,filenames.trimmedr,0.6,"Trimmed reads already filtered by Kraken2! Skipping...","") is False:
            kraken2(input_args,filenames,kraken2_path)
            if lib.file_exists_bool(filenames.k_unclassifiedf,filenames.trimmedf,0.6,"Forward trimmed reads successfully filtered by Kraken2! Skipping...","Taxonomic filtering with Kraken2 failed. Check the logs. Continuing with non-filtered reads.") is False \
            and lib.file_exists_bool(filenames.k_unclassifiedr,filenames.trimmedr,0.6,"Reverse trimmed reads successfully filtered by Kraken2! Skipping...","Taxonomic filtering with Kraken2 failed. Check the logs. Continuing with non-filtered reads.") is False:
                filenames.trimmedf = filenames.k_unclassifiedf
                filenames.trimmedr = filenames.k_unclassifiedr

    else:
        lib.print_h("Skipping taxonomic filtering of trimmed reads")    

    assembly_text = "   __________\n___/ Assembly \_________________________________________________________________"
    lib.print_t(assembly_text)

    assembly_path = os.path.join("assembly",input_args.array)
    lib.make_path(assembly_path)
    filenames.assembly_fasta = os.path.join(assembly_path,"scaffolds.fasta")

    # Checks whether to perform a hybrid assembly (with Nanopore flag `-n`) or isolate assembly
    if input_args.nanopore is not None:
        if lib.file_exists(filenames.assembly_fasta,"Reads already assembled! Skipping...","Performing hybrid assembly with Flye please...") is False:
            print("hello there")
            print("test0")
            filenames.nanopore_sam = os.path.join(assembly_path,filenames.nanopore.split(".")[0]+".sam")
            filenames.racon_assembly = os.path.join("racon","*.fa*") # find out actual output filename
            filenames.medaka_assembly = os.path.join("medaka","*.fa*") # find out actual output filename
            filenames.assembly_bam = os.path.join(assembly_path,filenames.assembly_fasta.split(".")[0]+".bam")
            filenames.assembly_sorted_bam = os.path.join(assembly_path,filenames.assembly_fasta.split(".")[0]+"_sorted.bam") 
            print("test1")
            assembly_hybrid(input_args,filenames,assembly_path)
            hybrid_polish(input_args,filenames,assembly_path)
    else:
        if lib.file_exists(filenames.assembly_fasta,"Reads already assembled! Skipping...","Performing short read assembly with SPADes") is False:
            assembly_short(input_args,filenames,assembly_path)

    lib.file_exists_exit(filenames.assembly_fasta,"Reads successfully assembled!","Assembly failed... check the logs and your inputs")
