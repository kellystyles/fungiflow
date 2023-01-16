import os
import sys
import subprocess
import datetime
import library as lib

"""
The assembly module of the Fungiflow pipeline.
"""

def find_adapters(input_args, filenames, trimmed_path):
    """
    Finds adapter sequences in raw read files.  
    
    Inputs:     forward and reverse raw Illumina paired-end FASTQ reads.
    Outputs:    FASTA file containing identified adapter sequences.
    """

    stdout = os.path.join(trimmed_path, f"adapter_{input_args.array}.out")
    stderr = os.path.join(trimmed_path, f"adapter_{input_args.array}.err")

    adapter_sh = os.path.join(os.path.dirname(os.path.realpath(__file__)), "adap_ID.sh")

    lib.print_h(f"Identifying adapter sequences in {input_args.array} reads")
    cmd = ["bash", adapter_sh, filenames.shortf, filenames.adapter_file]
    cmd2 = ["bash", adapter_sh, filenames.shortr, filenames.adapter_file]

    lib.execute(cmd, stdout, stderr)
    lib.execute(cmd2, stdout, stderr)

def trimmomatic(input_args, filenames, trimmed_path):
    """
    Trims primer sequences from raw Illumina reads with `trimmomatic`.  
    
    Inputs:     forward and reverse Illumina paired-end FASTQ reads.
    Outputs:    forward and reverse trimmed  Illumina paired-end FASTQ reads.
    """

    stdout = os.path.join(trimmed_path, f"{input_args.array}.out")
    stderr = os.path.join(trimmed_path, f"{input_args.array}.err")

    lib.print_h(f"Trimming Illumina {input_args.array} reads with Trimmomatic")
    cmd = ["TrimmomaticPE", "-threads", input_args.cpus, "-trimlog", filenames.trimlog_file, \
        filenames.shortf, filenames.shortr, filenames.trimmedf, filenames.utrimmedf, \
            filenames.trimmedr, filenames.utrimmedr, \
                #f"ILLUMINACLIP:{filenames.adapter_file}:2:30:10:4:4:/true", \
                    "TRAILING:9", "SLIDINGWINDOW:4:15", "MINLEN:36"]
    if len(filenames.singularity) > 0: cmd = filenames.singularity + cmd

    lib.execute(cmd, stdout, stderr)

def fastqc(input_args, filenames, trimmed_path):
    """
    Performs QC analysis of Illumina reads with `fastqc`.  
    
    Inputs:     forward and reverse Illumina paired-end FASTQ reads.
    Outputs:    fastqc output folder
    """

    stdout = os.path.join(trimmed_path, f"{input_args.array}_fastqc.out")
    stderr = os.path.join(trimmed_path, f"{input_args.array}_fastqc.err")
    lib.print_n("Assessing reads with FastQC")
    cmd = ["fastqc", "-o", trimmed_path, \
        input_args.illumina_f, input_args.illumina_r, \
        filenames.trimmedf, filenames.trimmedr]
    if len(filenames.singularity) > 0: cmd = filenames.singularity + cmd

    lib.execute(cmd, stdout, stderr)
    lib.file_exists_exit(filenames.fastqc_out, \
            "fastqc completed", "fastqc failed")

def porechop(input_args, filenames, trimmed_path):
    stdout = os.path.join(trimmed_path, f"{input_args.array}_porechop.out")
    stderr = os.path.join(trimmed_path, f"{input_args.array}_porechop.err")
    # Trimming long reads with PORECHOP
    cmd = ["porechop", "-i", filenames.nanopore, "-o", filenames.nanopore_trimmed, "--adapter_threshold", "96", "--no_split"]
    if len(filenames.singularity) > 0: cmd = filenames.singularity + cmd

    lib.execute(cmd, stdout, stderr)

def fastq_2_fasta(input_args, filenames, trimmed_path):
    stdout = os.path.join(trimmed_path, f"{input_args.array}_fastq_2_fasta.out")
    stderr = os.path.join(trimmed_path, f"{input_args.array}_fastq_2_fasta.err")
    # Converting long reads to FASTA format
    cmd = ["reformat.sh", f"in={filenames.nanopore_trimmed}", f"out={filenames.nanopore_fasta}", "ignorebadquality=True"] 
    if len(filenames.singularity) > 0: cmd = filenames.singularity + cmd

    lib.execute(cmd, stdout, stderr)

def length_filter_fasta(input_args, filenames, trimmed_path):
    stdout = os.path.join(trimmed_path, f"{input_args.array}_length_filter_fasta.out")
    stderr = os.path.join(trimmed_path, f"{input_args.array}_length_filter_fasta.err")
    # Length-filtering trimmed long reads
    cmd = ["bbduk.sh", f"in={filenames.nanopore_fasta}", f"out={filenames.nanopore_length_filtered}", f"minlen={input_args.minimum_length}"]
    if len(filenames.singularity) > 0: cmd = filenames.singularity + cmd

    lib.execute(cmd, stdout, stderr)

def short_map_2_long(input_args, filenames, trimmed_path):
    stdout = os.path.join(trimmed_path, f"{input_args.array}_short_map_2_long.out")
    stderr = os.path.join(trimmed_path, f"{input_args.array}_short_map_2_long.err")
    # mapping short reads
    a = f"gunzip -c {filenames.trimmedf} {filenames.trimmedr}"
    b = "awk \'NR % 4 == 2\'"
    c = "sort"
    d = "tr NT TN"
    e = "ropebwt2 -LR"
    if len(filenames.singularity) > 0: e = " ".join(filenames.singularity) + " " + e
    f = "tr NT TN"
    g = f"fmlrc-convert {filenames.short_mapped_2_long}"
    if len(filenames.singularity) > 0: g = " ".join(filenames.singularity) + " " + g
    commands = [a, b, c, d, e, f, g]

    lib.execute_shell(commands, stdout, stderr)

def long_read_corr(input_args, filenames, trimmed_path):
    stdout = os.path.join(trimmed_path, f"{input_args.array}_long_read_corr.out")
    stderr = os.path.join(trimmed_path, f"{input_args.array}_long_read_corr.err")
    # Correcting long reads with FMLRC
    cmd = ["fmlrc", "-p", input_args.cpus, filenames.short_mapped_2_long, filenames.nanopore_fasta, filenames.nanopore_corr]
    if len(filenames.singularity) > 0: cmd = filenames.singularity + cmd

    lib.execute(cmd, stdout, stderr)

def long_read_trim(input_args, filenames, trimmed_path):
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
    stdout = os.path.join(trimmed_path, f"{input_args.array}_long_read_trim.out")
    stderr = os.path.join(trimmed_path, f"{input_args.array}_long_read_trim.err")

    lib.print_h(f"Preparing long reads {input_args.nanopore} for assembly")
    if lib.file_exists(filenames.nanopore_trimmed, \
        "porechop already trimmed long reads. Skipping...", "Trimming long reads with porechop") == False:
        porechop(input_args, filenames, trimmed_path)
        lib.file_exists_exit(filenames.nanopore_trimmed, \
            "Porechop successfully trimmed the reads", "Porechop did not complete successfully. Check the logs.")
    if  lib.file_exists(filenames.nanopore_fasta, \
        "reformat.sh already converted the reads to FASTA format. Skipping...", "Using reformat.sh to convert long reads to FASTA format") == False:
        fastq_2_fasta(input_args, filenames, trimmed_path)
        lib.file_exists_exit(filenames.nanopore_fasta, \
            "reformat.sh successfully converted the reads to FASTA format", "reformat.sh did not complete successfully. Check the logs.")   
    if lib.file_exists(filenames.nanopore_length_filtered, \
        "bbduk.sh already length-filtered the reads. Skipping...", "Using bbduk.sh to length-filter the long reads") == False:
        length_filter_fasta(input_args, filenames, trimmed_path)
        lib.file_exists_exit(filenames.nanopore_length_filtered, \
            "bbduk.sh successfully length-filtered the reads", "bbduk.sh did not complete successfully. Check the logs.")
    if lib.file_exists(filenames.short_mapped_2_long, \
        "short reads already mapped to the long reads. Skipping...", "Using fmlrc and ropebwt2 to map the short reads to the long reads") == False:
        short_map_2_long(input_args, filenames, trimmed_path)
        lib.file_exists_exit(filenames.short_mapped_2_long, \
            "The short reads were successfully mapped to the long reads", "The short read mapping to long reads was not successful. Check the logs.")
    if lib.file_exists(filenames.nanopore_corr, \
        "fmlrc already corrected the reads", "Using fmlrc to correct the long reads") == False:
        long_read_corr(input_args, filenames, trimmed_path)
        lib.file_exists_exit(filenames.nanopore_corr, \
            "fmlrc successfully corrected the long reads", "fmlrc did not complete successfully. Check the logs.")
    lib.file_exists_exit(filenames.nanopore_corr, \
        "Long reads trimmed and corrected!", "Trimming and correction of long reads failed. Please check the logs")

def assembly_short(input_args, filenames, assembly_path):
    """
    Performs assembly of short reads with `spades` assembler.  
    
    Inputs:     forward and reverse trimmed Illumina paired-end FASTQ reads.
    Outputs:    spades output folder
    """

    stdout = os.path.join(assembly_path, f"{input_args.array}.out")
    stderr = os.path.join(assembly_path, f"{input_args.array}.err")

    # Defines the command based on input arguments
    if input_args.type == "metagenomic":
        cmd = ["spades.py", "--meta", "--threads", input_args.cpus, "--memory", input_args.mem, \
                "-k", "21, 33, 55, 77, 99, 127", \
                    "--pe1-1", filenames.trimmedf, "--pe1-2", filenames.trimmedr, "-o", assembly_path] 
    else:
        if input_args.careful is True:
            cmd = ["spades.py", "--threads", input_args.cpus, "--memory", input_args.mem, \
                     "--sc", "--pe1-1", filenames.trimmedf, "--pe1-2", filenames.trimmedr, "-o", assembly_path] 
        else:
            cmd = ["spades.py", "--isolate", "--threads", input_args.cpus, "--memory", input_args.mem, \
                    "-k", "21, 33, 55, 77, 99, 127", "--pe1-1", filenames.trimmedf, \
                        "--pe1-2", filenames.trimmedr, "-o", assembly_path] 
    # Assembly with SPADES
    if len(filenames.singularity) > 0: cmd = filenames.singularity + cmd
    lib.print_h(f"Assembling trimmed {input_args.array} reads with SPADES")

    lib.execute(cmd, stdout, stderr)
    lib.file_exists_exit(filenames.assembly_fasta, \
        "Short reads successfully assembled with SPADes", "SPADes assembly failed... check the logs and your inputs")


def assembly_hybrid(input_args, filenames, assembly_path):
    """
    Performs a hybrid assembly using `flye` assembler.

    Inputs:     forward and reverse trimmed Illumina paired-end FASTQ reads, trimmed Nanopore FASTA reads.
    Outputs:    flye output folder
    """
    stdout = os.path.join(assembly_path, f"{input_args.array}.out")
    stderr = os.path.join(assembly_path, f"{input_args.array}.err")

    if input_args.type == "metagenomic":
        cmd = ["flye", "--nano-raw", filenames.nanopore_corr, "--out-dir", assembly_path, "--threads", str(input_args.cpus), "--iterations", "2", "--meta"]
    else:
        cmd = ["flye", "--nano-raw", filenames.nanopore_corr, "--out-dir", assembly_path, "--threads", str(input_args.cpus), "--iterations", "2"]
    if len(filenames.singularity) > 0: cmd = filenames.singularity + cmd    
    # Assembly with FLYE
    lib.print_h(f"Hybrid isolate assembly of {input_args.array} reads with FLYE")

    lib.execute(cmd, stdout, stderr)
    lib.file_exists_exit(filenames.assembly_fasta, \
        "Reads successfully assembled with flye!", "Assembly with flye failed... check the logs and your inputs")
    print(f"{datetime.datetime.now():%Y-%m-%d %I:%M:%S} Completed assembly of trimmed {input_args.array} reads with FLYE")

def minimap2(input_args, filenames, assembly_path, read_length):
    stdout = os.path.join(assembly_path, f"{input_args.array}_minimap2.out")
    stderr = os.path.join(assembly_path, f"{input_args.array}_minimap2.err")
    
    if read_length == "long":
        cmd = ["minimap2", "-ax", "map-ont", filenames.assembly_fasta, filenames.nanopore_corr, "-o", filenames.nanopore_sam]
    elif read_length == "short":
        cmd = ["minimap2", "-ax", "sr", filenames.medaka_consensus, filenames.trimmedf, filenames.trimmedr, "-o", filenames.medaka_consensus_sam]
    else:
        print("\'read_length\' parameter unknown. Accepted values are \"short\" or \"long\"")
        sys.exit()
    if len(filenames.singularity) > 0: cmd = filenames.singularity + cmd    

    lib.execute(cmd, stdout, stderr)
    if read_length == "long":
        lib.file_exists_exit(filenames.nanopore_sam, \
            f"Mapping of {filenames.nanopore_length_filtered} to {filenames.assembly_fasta} with minimap2 was successful", f"Mapping of {filenames.nanopore_length_filtered} to {filenames.assembly_fasta} with minimap2 failed.")
    if read_length == "short":
        lib.file_exists_exit(filenames.medaka_consensus_sam, \
            f"Mapping of {filenames.trimmedf} & {filenames.trimmedr} to {filenames.medaka_consensus} with minimap2 was successful", f"Mapping of {filenames.trimmedf} & {filenames.trimmedr} to {filenames.medaka_consensus} with minimap2 failed")

def racon(input_args, filenames, assembly_path):
    stdout = filenames.racon_consensus
    stderr = os.path.join(assembly_path, f"{input_args.array}_racon.out")
    cmd = ["racon", "--no-trimming", "-t", input_args.cpus, filenames.nanopore_corr, filenames.nanopore_sam, filenames.assembly_fasta]
    if len(filenames.singularity) > 0: cmd = filenames.singularity + cmd
   
    lib.execute(cmd, stdout, stderr)
    lib.file_exists_exit(filenames.racon_consensus, \
        f"Polishing {filenames.assembly_fasta} with racon successful", f"Polishing {filenames.assembly_fasta} with racon failed")

def medaka(input_args, filenames, assembly_path):
    stdout = os.path.join(assembly_path, f"{input_args.array}_medaka.out")
    stderr = os.path.join(assembly_path, f"{input_args.array}_medaka.err")
    cmd1 = ["mini_align", "-i", filenames.nanopore_length_filtered, "-r", filenames.racon_consensus, "-m", "-p", filenames.racon_consensus, "-t", input_args.cpus]
    cmd2 = ["medaka", "consensus", filenames.racon_consensus_bam, filenames.racon_consensus_hdf, "--threads", input_args.cpus]
    cmd3 = ["medaka", "stitch", filenames.racon_consensus_hdf, filenames.racon_consensus, filenames.medaka_consensus]
    if len(filenames.singularity) > 0: 
        cmd1 = filenames.singularity + cmd1
        cmd2 = filenames.singularity + cmd2
        cmd3 = filenames.singularity + cmd3

    lib.execute(cmd1, stdout, stderr)
    lib.file_exists_exit(filenames.racon_consensus_bam, \
        "mini_align completed successfully", "mini_align failed")
    lib.execute(cmd2, stdout, stderr)
    lib.file_exists_exit(filenames.racon_consensus_hdf, \
        "medaka consensus completed successfully", "medaka consensus failed")
    lib.execute(cmd3, stdout, stderr)
    lib.file_exists_exit(filenames.medaka_consensus, \
        f"medaka stitch completed successfully\nPolishing {filenames.racon_consensus} with medaka successful", "medaka stitch failed")

def samtools_sort(input_args, filenames, assembly_path):
    stdout = os.path.join(assembly_path, f"{input_args.array}_samtools_sort.out")
    stderr = os.path.join(assembly_path, f"{input_args.array}_samtools_sort.err")
    cmd = ["samtools", "sort", "-@", input_args.cpus, "-o", filenames.medaka_sorted_sam, filenames.medaka_consensus_sam]
    if len(filenames.singularity) > 0: cmd = filenames.singularity + cmd

    lib.execute(cmd, stdout, stderr)
    lib.file_exists(filenames.medaka_sorted_sam, \
        f"{filenames.medaka_consensus_sam} successfully sorted with samtools sort", f"Sorting {filenames.medaka_consensus_sam} with samtools failed")

def samtools_index(input_args, filenames, assembly_path):
    stdout = os.path.join(assembly_path, f"{input_args.array}_samtools_index.out")
    stderr = os.path.join(assembly_path, f"{input_args.array}_samtools_index.err")
    cmd = ["samtools", "index", "-@", input_args.cpus, filenames.medaka_sorted_sam]
    if len(filenames.singularity) > 0: cmd = filenames.singularity + cmd

    lib.execute(cmd, stdout, stderr)
    lib.file_exists(filenames.assembly_sorted_index, \
        f"{filenames.medaka_sorted_sam} successfully indexed with samtools index", f"Indexing {filenames.medaka_sorted_sam} indexing with samtools failed")

def polypolish(input_args, filenames, assembly_path):
    """
    Polishes an assembly with polypolish.
    Input: medaka consensus FASTA, short reads mapped and sorted to medaka consensus SAM
    Output: polypolish consensus FASTA 
    """
    stdout = os.path.join(assembly_path, f"{input_args.array}_polypolish.out")
    stderr = os.path.join(assembly_path, f"{input_args.array}_polypolish.err")
    cmd = ["polypolish", filenames.medaka_consensus, filenames.medaka_sorted_sam, ">", filenames.polypolish_consensus]
    if len(filenames.singularity) > 0: cmd = filenames.singularity + cmd

    lib.execute(cmd, stdout, stderr)
    lib.file_exists(filenames.polypolish_consensus, \
        f"Polishing {filenames.medaka_consensus} with polypolish was successful", f"Polishing {filenames.medaka_consensus} with polypolish failed")    

def hybrid_polish(input_args, filenames, assembly_path):
    """
    Polishes hybrid assemblies as follows:
    1. map trimmed/length-filtered long reads to assembly with `minimap2`
    2. one round of polishing using `racon`
    3. one round of polishing using `medaka`
    4. map short trimmed/taxa-filtered short reads to medaka consensus assembly
    5. one round of polishing using `polypolish`
    
    Inputs:    hybrid flye assembly FASTA file, Nanopore FASTA and Illumina FASTQ reads.
    Output:    polished assembly FASTA file.
    """
    stdout = os.path.join(assembly_path, f"{input_args.array}_polish.out")
    stderr = os.path.join(assembly_path, f"{input_args.array}_polish.err")
    
    # Polishing FLYE assembly
    if lib.file_exists(filenames.nanopore_sam, \
        "Nanopore SAM file already exists, Skipping...", f"Mapping {filenames.nanopore_length_filtered} reads to {filenames.assembly_fasta} with minimap2") == False:
        minimap2(input_args, filenames, assembly_path, "long")
    if lib.file_exists(filenames.racon_consensus, \
        "racon consensus assembly already exists, Skipping...", f"Polishing {filenames.assembly_fasta} with racon") == False:
        racon(input_args, filenames, assembly_path)
    if lib.file_exists(filenames.medaka_consensus, \
        "medaka consensus assembly already exists, Skipping...", f"Polishing {filenames.racon_consensus} with medaka") == False:
        medaka(input_args, filenames, assembly_path)
    if lib.file_exists(filenames.medaka_consensus_sam, \
        f"{filenames.medaka_consensus_sam} already exists. Skipping...", f"Mapping {filenames.trimmedf} & {filenames.trimmedr} to {filenames.medaka_consensus} with minimap2...") == False:
        minimap2(input_args, filenames, assembly_path, "short")
    if lib.file_exists(filenames.medaka_sorted_sam, \
        f"{filenames.medaka_sorted_sam} already exists. Skipping...", f"Sorting {filenames.medaka_consensus_sam} with samtools...") == False:
        samtools_sort(input_args, filenames, assembly_path)
    if lib.file_exists(filenames.assembly_sorted_index, \
        f"{filenames.assembly_sorted_index} already exists. Skipping...", f"Indexing {filenames.medaka_sorted_sam} with samtools...") == False:
        samtools_index(input_args, filenames, assembly_path)
    if lib.file_exists(filenames.polypolish_consensus, \
        f"{filenames.polypolish_consensus} consensus already exists. Skipping...", f"Polishing {filenames.medaka_sorted_sam} with polypolish...") == False:
        polypolish(input_args, filenames, assembly_path)

def kraken2(input_args, filenames, kraken_path, read_length):
    """
    Taxonomically filters trimmed Illumina reads using `kraken2`.

    Inputs:     forward and reverse trimmed Illumina reads.
    Outputs:    forward and reverse filtered and unfiltered Illumina reads.
    """
    
    stdout = os.path.join(kraken_path, f"{input_args.array}.out")
    stderr = os.path.join(kraken_path, f"{input_args.array}.err")

    if read_length == "short":
        lib.print_h(f"Performing Kraken2 analysis of short trimmed reads {input_args.array}")
        cmd = ["kraken2", "--db", input_args.kraken2_db, \
            "--threads", input_args.cpus, "--use-names", "--report", filenames.kreport, "--classified-out", filenames.k_classified, \
                "--unclassified-out", filenames.k_unclassified, "--paired", filenames.trimmedf, filenames.trimmedr]
    elif read_length == "long":
        lib.print_h(f"Performing Kraken2 analysis of long trimmed reads {input_args.array}")
        cmd = ["kraken2", "--db", input_args.kraken2_db, \
            "--threads", input_args.cpus, "--use-names", "--report", filenames.kreport_long, "--classified-out", filenames.k_long_classified, \
                "--unclassified-out", filenames.k_long_unclassified, filenames.nanopore_corr]
    else:
        lib.print_e("The `read_length` parameter for the kraken2 function is not valid.")
    if len(filenames.singularity) > 0: cmd = filenames.singularity + cmd
        
    lib.execute(cmd, stdout, stderr)

def main(input_args, filenames):

    lib.print_h("Initializing \'assembly\' module...")
    assembly_start_time = datetime.datetime.now()
    trim_text = "  __________\n___/ Trimming \_________________________________________________________________"
    lib.print_t(trim_text)

    trimmed_path = os.path.join(input_args.directory_new, "trimmed")
    lib.make_path(trimmed_path)

    #?# if directory input_arg == '.' then this command will put the input_arg.array value before the read, so need to use absolute path for the directory argument
    #filenames.adapter_file = os.path.join(trimmed_path, f"{str(input_args.array)}_adapters.fasta")
    #if lib.file_exists(filenames.adapter_file, "Adapters already identified! Skipping...", "Searching for adapter sequences in short reads") is False:
     #   find_adapters(input_args, filenames, trimmed_path)
      #  lib.file_exists_exit(filenames.adapter_file, "Adapters identified!", "Adapters not successfully extracted... check the logs and your inputs")

    filenames.trimmedf = os.path.join(trimmed_path, f"{input_args.array}_trimmed_1P.fq.gz")
    filenames.trimmedr = os.path.join(trimmed_path, f"{input_args.array}_trimmed_2P.fq.gz")
    filenames.utrimmedf = os.path.join(trimmed_path, f"{input_args.array}_trimmed_1U.fq.gz")
    filenames.utrimmedr = os.path.join(trimmed_path, f"{input_args.array}_trimmed_2U.fq.gz")
    filenames.trimlog_file = os.path.join(trimmed_path, f"{input_args.array}_trimmomatic.log")

    if lib.file_exists_bool(filenames.trimmedf, filenames.shortf, 0.6, "Forward reads already trimmed! Skipping...", "") == False and \
        lib.file_exists_bool(filenames.trimmedr, filenames.shortr, 0.6, "Reverse reads already trimmed! Skipping...", "") == False:
        trimmomatic(input_args, filenames, trimmed_path)
        lib.file_exists_exit(filenames.trimmedf, "Forward reads trimmed successfully!", "Forward reads not trimmed... check the logs and your inputs")
        lib.file_exists_exit(filenames.trimmedr, "Reverse reads trimmed successfully!", "Reverse reads not trimmed... check the logs and your inputs")

    if input_args.nanopore is not None:
        filenames.nanopore = os.path.abspath(os.path.join(input_args.directory, input_args.nanopore))
        filenames.nanopore_trimmed = os.path.join(trimmed_path, f"{input_args.array}_ONT_trimmed.fq")
        filenames.nanopore_fasta = os.path.join(trimmed_path, f"{input_args.array}_ONT_trimmed.fa")
        filenames.nanopore_length_filtered = os.path.join(trimmed_path, f"{input_args.array}_lf.fa")
        filenames.short_mapped_2_long = os.path.join(trimmed_path, f"{input_args.array}_mapped.npy")
        filenames.nanopore_corr = os.path.join(trimmed_path, f"{input_args.array}_corr.fa")
        if lib.file_exists(filenames.nanopore_corr, "Long reads already trimmed and corrected! Skipping...", "") == False:
            long_read_trim(input_args, filenames, trimmed_path)

    #filenames.fastqc_out = os.path.join(trimmed_path, filenames.trimmedf.split(".")[0]+"_fastqc.zip")
    #if lib.file_exists(filenames.fastqc_out, "", "") == False:
     #   fastqc(input_args, filenames, trimmed_path)

    if input_args.kraken2 is True:

        tax_filter_text = "  _____________________\n___/ Taxonomic Filtering \______________________________________________________"
        lib.print_t(tax_filter_text)

        kraken2_path = os.path.join(input_args.directory_new, "kraken2")
        lib.make_path(kraken2_path)

        # short reads
        filenames.kreport = os.path.join(kraken2_path, f"{input_args.array}_short.report")
        filenames.k_classified = os.path.join(kraken2_path, f"{input_args.array}_class#.fq")
        filenames.k_unclassified = os.path.join(kraken2_path, f"{input_args.array}_unclass#.fq")
        filenames.k_classifiedf = os.path.join(kraken2_path, f"{input_args.array}_class_1.fq")
        filenames.k_classifiedr = os.path.join(kraken2_path, f"{input_args.array}_class_2.fq")
        filenames.k_unclassifiedf = os.path.join(kraken2_path, f"{input_args.array}_unclass_1.fq")
        filenames.k_unclassifiedr = os.path.join(kraken2_path, f"{input_args.array}_unclass_2.fq")

        if lib.file_exists_bool(filenames.k_unclassifiedf, filenames.trimmedf, 0.05, "Short trimmed reads already filtered by Kraken2! Skipping...", "") == False \
        and lib.file_exists_bool(filenames.k_unclassifiedr, filenames.trimmedr, 0.05, "Short trimmed reads already filtered by Kraken2! Skipping...", "") == False:
            kraken2(input_args, filenames, kraken2_path, "short")
            if lib.file_exists_bool(filenames.k_unclassifiedf, filenames.trimmedf, 0.05, "Forward short trimmed reads successfully filtered by Kraken2!", "Taxonomic filtering short reads with Kraken2 failed. Check the logs. Continuing with non-filtered reads.") == True \
            and lib.file_exists_bool(filenames.k_unclassifiedr, filenames.trimmedr, 0.05, "Reverse short trimmed reads successfully filtered by Kraken2!", "Taxonomic filtering short reads with Kraken2 failed. Check the logs. Continuing with non-filtered reads.") == True:
                filenames.trimmedf = filenames.k_unclassifiedf
                filenames.trimmedr = filenames.k_unclassifiedr
        else:
            filenames.trimmedf = filenames.k_unclassifiedf
            filenames.trimmedr = filenames.k_unclassifiedr

    else:
        lib.print_h("Skipping taxonomic filtering of trimmed reads")    

    assembly_text = "  __________\n___/ Assembly \_________________________________________________________________"
    lib.print_t(assembly_text)

    assembly_path = os.path.join(input_args.directory_new, "assembly")
    lib.make_path(assembly_path)

    # Checks whether to perform a hybrid assembly (with Nanopore flag `-n`) or isolate assembly
    if input_args.nanopore is not None:
        filenames.assembly_fasta = os.path.join(assembly_path, "assembly.fasta")
        filenames.nanopore_sam = os.path.join(assembly_path, f"{input_args.array}_ONT_corr.sam")
        filenames.racon_consensus = os.path.join(assembly_path, "racon_consensus.fa")
        filenames.racon_consensus_bam = os.path.join(assembly_path, "racon_consensus.fa.bam")
        filenames.racon_consensus_index = os.path.join(assembly_path, "racon_consensus.fa.bam.bai")
        filenames.racon_consensus_hdf = os.path.join(assembly_path, "racon_consensus.hdf")
        filenames.medaka_consensus = os.path.join(assembly_path, "medaka_consensus.fa")
        filenames.medaka_index = os.path.join(assembly_path, "medaka_consensus.bai")
        filenames.medaka_consensus_sam = os.path.join(assembly_path, "medaka_consensus.sam")
        filenames.medaka_sorted_sam = os.path.join(assembly_path, "medaka_sorted.sam")
        filenames.assembly_sorted_index = os.path.join(assembly_path, "medaka_sorted.sam.bai")
        filenames.polypolish_consensus = os.path.join(assembly_path, "polypolish.fasta")
        if lib.file_exists(filenames.assembly_fasta, "Reads already assembled! Skipping...", "Performing hybrid assembly with Flye...") == False: 
            assembly_hybrid(input_args, filenames, assembly_path)
        polishing_text = "  ___________\n___/ Polishing \_________________________________________________________________"
        lib.print_t(polishing_text)
        if lib.file_exists(filenames.polypolish_consensus, "Hybrid assembly already polished! Skipping...", "Polishing hybrid assembly...") == False:
            hybrid_polish(input_args, filenames, assembly_path)
            filenames.assembly_fasta = str(filenames.polypolish_consensus)
    else:
        filenames.assembly_fasta = os.path.join(assembly_path, "scaffolds.fasta")
        if lib.file_exists(filenames.assembly_fasta, "Reads already assembled! Skipping...", "Performing short read assembly with SPADes") == False:
            assembly_short(input_args, filenames, assembly_path)

    lib.print_h(f"assembly module completed in {datetime.datetime.now() - assembly_start_time}")

if __name__ == '__main__':
    main()
