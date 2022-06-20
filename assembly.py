import os
import subprocess
import library as lib

"""
The assembly module of the Fungiflow pipeline.
"""

def find_adapters(input_args,filenames):
    """
    Finds adapter sequences in raw read files.  
    
    Inputs:     forward and reverse raw Illumina paired-end FASTQ reads.
    Outputs:    FASTA file containing identified adapter sequences.
    """

    stdout = os.path.join(filenames.adapter_path,f"adapter_{input_args.array}.out")
    stderr = os.path.join(filenames.adapter_path,f"adapter_{input_args.array}.err")

    print_h(f"Identifying adapter sequences in {input_args.array} reads")
    cmd1 = ["singularity","exec","-B","/nfs:/nfs",input_args.singularity,"bash","adap_ID.sh",filenames.forward_reads,filenames.adapter_file]
    cmd2 = ["singularity","exec","-B","/nfs:/nfs",input_args.singularity,"bash","adap_ID.sh",filenames.reverse_reads,filenames.adapter_file]
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

    print_h(f"Trimming Illumina {input_args.array} reads with Trimmomatic")
    cmd1 = ["singularity","exec","-B","/nfs:/nfs",input_args.singularity,"TrimmomaticPE","-threads",filenames.cpus,"-trimlog",filenames.log_file,filenames.forward_reads,filenames.reverse_reads,filenames.trimmedf,filenames.utrimmedf,filenames.trimmedr,filenames.utrimmedr,f"ILLUMINACLIP:{filenames.adapter_file}:2:30:10:4:4:/true","TRAILING:9","SLIDINGWINDOW:4:15","MINLEN:36"]
    print(" ".join(cmd1))
    try:
        lib.execute(cmd1,stdout,stderr)
    except subprocess.CalledProcessError as e:
        print(e.returncode)
        print(e.output)

def fastqc(input_args,fastqc_path):
    """
    Performs QC analysis of Illumina reads with `fastqc`.  
    
    Inputs:     forward and reverse Illumina paired-end FASTQ reads.
    Outputs:    fastqc output folder
    """

    stdout = os.path.join(fastqc_path,f"{array}.out")
    stderr = os.path.join(fastqc_path,f"{array}.err")

    print_n("Assessing reads with FastQC")
    cmd1 = ["singularity","exec","-B","/nfs:/nfs",input_args.qc_image,"fastqc","-o",fastqc_path,f"{input_args.array}_*.fq*"]
    print(" ".join(cmd1))
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

    print_h(f"Preparing long reads, {input_args.nanopore}, for assembly")
    # Trimming long reads with PORECHOP
    filenames.nanopore_trimmed = os.path.join("trimmed",filenames.nanopore.split(".")[0]+".fq")
    cmd1 = ["singularity","exec","-B","/nfs:/nfs",input_args.singularity,"porechop","-i",filenames.nanopore,"-o",filenames.nanopore_trimmed,"--adapter_threshold","96","--no_split"]
    # Length-filtering trimmed long reads
    filenames.nanopore_length_filtered = os.path.join("trimmed",filenames.nanopore.split(".")[0]+"_lf.fq")
    cmd2 = ["singularity","exec","-B","/nfs:/nfs",input_args.singularity,"bbduk.sh",f"in={input_args.nanopore_trimmed}",f"out={filenames.nanopore_length_filtered}","minlen=3000"]
    # Converting long reads to FASTA format
    filenames.nanopore_fasta = os.path.join("trimmed",filenames.nanopore.split(".")[0]+"_lf.fasta")
    cmd3 = ["singularity","exec","-B","/nfs:/nfs",input_args.singularity,"reformat.sh",f"in={input_args.nanopore_length_filtered}",f"out={filenames.nanopore_fasta}"]
    # mapping short reads
    filenames.short_mapped_2_long = os.path.join("trimmed",input_args.array+"mapped.npy")
    cmd4 = ["gunzip","-c",filenames.trimmedf,filenames.trimmedr,"|","awk","\'NR % 4 == 2\'","|","sort","|","tr","NT","TN","|","singularity","exec","-B","/nfs:/nfs",input_args.singularity,"ropebwt2","-LR","|","tr","NT","TN","|","singularity","exec","-B","/nfs:/nfs",input_args.singularity,"fmlrc-convert",filenames.short_mapped_2_long]
    # Correcting long reads with FMLRC
    filenames.nanopore_corr = os.path.join("trimmed",filenames.nanopore.split(".")[0]+"corr.fasta")
    cmd5 = ["singularity","exec","-B","/nfs:/nfs",input_args.singularity,"fmlrc","-p",input_args.cpus,filenames.short_mapped_2_long,filenames.nanopore_fasta,filenames.nanopore_corr]
    
    print(" ".join(cmd1))
    print(" ".join(cmd2))
    print(" ".join(cmd3))
    print(" ".join(cmd4))
    print(" ".join(cmd5))
    try:
        lib.execute(cmd1,stdout,stderr)
        lib.execute(cmd2,stdout,stderr)
        lib.execute(cmd3,stdout,stderr)
        lib.execute_shell(cmd4,stdout,stderr)
        lib.execute(cmd5,stdout,stderr)
    except subprocess.CalledProcessError as e:
        print(e.returncode)
        print(e.output)

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
    print_h(f"Assembling trimmed {input_args.array} reads with SPADES")
    cmd1 = ["singularity","exec","-B","/nfs:/nfs",input_args.singularity,"spades.py",assembly_type,"--threads",input_args.cpus,"--memory",input_args.mem,"-k","21,33,55,77,99,127","--pe1-1",filenames.trimmedf,"--pe1-2",filenames.trimmedr,"-o",assembly_path]
    print(" ".join(cmd1))
    try:
        lib.execute(cmd1,stdout,stderr)
        print(f"{datetime.datetime.now():%Y-%m-%d %I:%M:%S} Completed assembly of trimmed {input_args.array} reads with SPADES")
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

    if input_args.data_type is "metagenomic":
        assembly_type = "--meta"
    else:
        assembly_type = ""

    # Assembly with FLYE
    print_h(f"Hybrid isolate assembly of {input_args.array} reads with FLYE")
    cmd1 = ["singularity","exec","-B","/nfs:/nfs",input_args.singularity,"flye","--nano-raw",filenames.nanopore_clean,"-o",assembly_path,"-t",input_args.cpus,"--trestle","-i","2",assembly_type]
    print(" ".join(cmd1))
    try:
        lib.execute(cmd1,stdout,stderr)
        print(f"{datetime.datetime.now():%Y-%m-%d %I:%M:%S} Completed assembly of trimmed {input_args.array} reads with SPADES")
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
    
    Inputs:     hybrid flye assembly FASTA file, Nanopore FASTA and Illumina FASTQ reads.
    Output:    polished assembly FASTA file.
    """

    # Polishing FLYE assembly
    print_h(f"Polishing {input_args.assembly_fasta}")
    filenames.nanopore_sam = os.path.join(assembly_path,filenames.nanopore.split(".")[0]+".sam")
    cmd1 = ["singularity","exec","-B","/nfs:/nfs",input_args.singularity,"minimap2","-x","map-ont",filenames.assembly_fasta,filenames.nanopore_length_filtered,">",filenames.nanopore_sam]
    racon_path = os.path.join(assembly_path,"racon")
    lib.make_path(racon_path)
    os.chdir(racon_path)
    cmd2 = ["singularity","exec","-B","/nfs:/nfs",input_args.singularity,"racon","--no-trimming","-t",input_args.cpus,filenames.nanopore_length_filtered,filenames.nanopore_sam,filenames.assembly_fasta]
    filenames.racon_assembly = os.path.join(racon_path,"*.fa*") # find out actual output filename
    medaka_path = os.path.join(assembly_path,"medaka")
    lib.make_path(medaka_path)
    os.chdir(medaka_path)
    filenames.medaka_assembly = os.path.join(medaka_path,"*.fa*") # find out actual output filename
    cmd3 = ["singularity","exec","-B","/nfs:/nfs",input_args.singularity,"medaka_consenus","-i",filenames.nanopore_length_filtered,"-d",filenames.racon_assembly,"-o","medaka","-t",input_args.cpus,"-m","r941_min_hac_g511"]
    os.chdir(input_args.directory)
    cmd4 = ["singularity","exec","-B","/nfs:/nfs",input_args.singularity,"bwa-mem2","index",filenames.assembly_fasta]
    filenames.assembly_bam = os.path.join(assembly_path,filenames.assembly_fasta.split(".")[0]+".bam")
    cmd4 = ["singularity","exec","-B","/nfs:/nfs",input_args.singularity,"bwa-mem2","mem","-t",input_args.cpus,filenames.medaka_assembly,filenames.trimmedf,filenames.trimmedr,">",filenames.assembly_bam]
    filenames.assembly_sorted_bam = os.path.join(assembly_path,filenames.assembly_fasta.split(".")[0]+"_sorted.bam")    
    cmd5 = ["singularity","exec","-B","/nfs:/nfs",input_args.singularity,"samtools","sort","-@",input_args.cpus,"-o",filenames.assembly_sorted_bam,filenames.assembly_bam]
    cmd6 = ["singularity","exec","-B","/nfs:/nfs",input_args.singularity,"samtools","index","-@",input_args.cpus,filenames.assembly_sorted_bam]
    cmd7 = ["singularity","exec","-B","/nfs:/nfs",input_args.singularity,"java",f"-Xmx{input_args.mem}G","-jar","pilon","--genome",input_args.assembly_fasta,"--bam",filenames.assembly_sorted_bam,"--outdir",assembly_path]

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

def kraken2(input_args,filenames):
    """
    Taxonomically filters trimmed Illumina reads using `kraken2`.

    Inputs:     forward and reverse trimmed Illumina reads.
    Outputs:    forward and reverse filtered and unfiltered Illumina reads.
    """
    
    stdout = os.path.join(kraken_path,f"{array}.out")
    stderr = os.path.join(kraken_path,f"{array}.err")

    print_h(f"Performing Kraken2 analysis of trimmed reads for {input_args.array}")
    cmd1 = ["singularity","exec","-B","/nfs:/nfs",input_args.singularity,"kraken2","--db",input_args.kraken2_db,\
        "--threads",input_args.cpus,"--use-filenames","--report",filenames.kreport,"--classified-out",classified,\
            "--unclassified-out",unclassified,"--paired",trimmed_forward,trimmed_reverse]
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
    if lib.file_exists_bool(filenames.adapter_file,0,"Adapters already identified! Skipping...",) is False:
        lib.find_adapters(input_args,filenames)
        lib.lib.file_exists(filenames.adapter_file,0,"Adapters identified!","Adapters not successfully extracted... check the logs and your inputs")

    trimmed_path = "trimmed"
    lib.make_path(trimmed_path)

    filenames.trimmedf = os.path.join(trimmed_path,f"{input_args.array}_forward_trimmed_1P.fq.gz")
    filenames.trimmedr = os.path.join(trimmed_path,f"{input_args.array}_reverse_trimmed_2P.fq.gz")
    filenames.utrimmedf = os.path.join(trimmed_path,f"{input_args.array}_forward_trimmed_1U.fq.gz")
    filenames.utrimmedr = os.path.join(trimmed_path,f"{input_args.array}_reverse_trimmed_2U.fq.gz")
    filenames.trimlog_file = os.path.join(trimmed_path,f"{input_args.array}_trimmomatic.log")

    if lib.file_exists_bool(filenames.trimmedf,os.stat(filenames.trimmedf).st_size > (0.6 * os.stat(input_args.illumina_f).st_size),"Forward reads already trimmed! Skipping...",) is False and \
        lib.file_exists_bool(filenames.trimmedr,os.stat(filenames.trimmedr).st_size > (0.6 * os.stat(input_args.illumina_r).st_size),"Reverse reads already trimmed! Skipping...",) is False:
        lib.trimmomatic(input_args,filenames,trimmed_path)
        lib.file_exists(filenames.trimmedf,os.stat(filenames.trimmedf).st_size > (0.6 * os.stat(input_args.illumina_f).st_size),"Forward reads trimmed successfully!","Forward reads not trimmed... check the logs and your inputs")
        lib.file_exists(filenames.trimmedr,os.stat(filenames.trimmedr).st_size > (0.6 * os.stat(input_args.illumina_r).st_size),"Reverse reads trimmed successfully!","Reverse reads not trimmed... check the logs and your inputs")

    if input_args.type is "hybrid":
        long_read_trim(input_args,filenames,trimmed_path)

    if len(input_args.kraken2_db) > 0:

        tax_filter_text = "   _____________________\n___/ Taxonomic Filtering \______________________________________________________"
        lib.print_t(tax_filter_text)

        kraken2_path = "kraken2"
        lib.make_path(kraken2_path)
        filenames.kreport = os.path.join(kraken2_path,f"{input_args.array}.report")
        filenames.k_classifiedf = os.path.join(kraken2_path,f"{input_args.array}_class_1.fq")
        filenames.k_classifiedr = os.path.join(kraken2_path,f"{input_args.array}_class_2.fq")
        filenames.k_unclassifiedf = os.path.join(kraken2_path,f"{input_args.array}_unclass_1.fq")
        filenames.k_unclassifiedr = os.path.join(kraken2_path,f"{input_args.array}_unclass_2.fq")

        if lib.file_exists_bool(kraken2_path,os.stat(filenames.k_unclassifiedf).st_size > (0.6 * os.stat(filenames.trimmedf).st_size),"Forward trimmed reads already filtered by Kraken2! Skipping...",) is False \
        and lib.file_exists_bool(kraken2_path,os.stat(filenames.k_unclassifiedr).st_size > (0.6 * os.stat(filenames.trimmedr).st_size),"Reverse trimmed reads already filtered by Kraken2! Skipping...",) is False:
            kraken2(input_args,filenames)
            if lib.file_exists_bool(kraken2_path,os.stat(filenames.k_unclassifiedf).st_size > (0.6 * os.stat(filenames.trimmedf).st_size),"Forward trimmed reads successfully filtered by Kraken2!","Forward trimmed reads not filtered... check the logs and your inputs") is True \
            and lib.file_exists_bool(kraken2_path,os.stat(filenames.k_unclassifiedr).st_size > (0.6 * os.stat(filenames.trimmedr).st_size),"Reverse trimmed reads successfully filtered by Kraken2!","Reverse trimmed reads not filtered... check the logs and your inputs") is True:
                filenames.trimmedf = filenames.k_unclassifiedf
                filenames.trimmedr = filenames.k_unclassifiedr

    else:
        lib.print_h("Skipping taxonomic filtering of trimmed reads")    

    assembly_text = "   __________\n___/ Assembly \_________________________________________________________________"
    lib.print_t(assembly_text)

    assembly_path = os.path.join("assembly",input_args.array)
    lib.make_path(assembly_path)
    filenames.assembly_fasta = os.path.join(assembly_path,"scaffolds.fasta")

    if input_args.type is "short":
        if lib.file_exists_bool(filenames.assembly_fasta,0,"Reads already assembled! Skipping...",) is False:
            assembly_short(input_args,filenames,assembly_path)
    elif input_args.type is "hybrid":
        if lib.file_exists_bool(filenames.assembly_fasta,0,"Reads already assembled! Skipping...",) is False:
            assembly_hybrid(input_args,filenames,assembly_path)
            hybrid_polish(input_args,filenames,assembly_path)

    lib.file_exists(filenames.assembly_fasta,0,"Reads successfully assembled!","Assembly failed... check the logs and your inputs")
