import os
import subprocess
import datetime
import shutil
import pandas as pd
import library as lib

def itsx(input_args,filenames,itsx_path):
    """
    Extracts the ITS sequences from an assembly FASTA file using `itsx`.
    
    Input:      assembly FASTA file
    Output:     itsx directory with ITS FASTA files
    """
    stdout = os.path.join(itsx_path,f"{input_args.array}.out")
    stderr = os.path.join(itsx_path,f"{input_args.array}.err")

    cmd1 = ["singularity","exec","-B","/nfs:/nfs",input_args.singularity,"ITSx","-i",filenames.assembly_fasta,"-o",os.path.join(itsx_path,input_args.array),"-t","F","--nhmmer","T","--cpu",input_args.cpus,"--preserve","--only_full","T","--temp",itsx_path]
    try:
        lib.print_n("Extracting the ITS sequence with ITSx")
        lib.execute(cmd1,stdout,stderr)
    except subprocess.CalledProcessError as e:
        print(e.returncode)
        print(e.output)  

def blastn(input_args,filenames,blast_path):
    """
    Will run `blastn` on an ITS FASTA file and output to a results TXT file.

    Input:      ITS input FASTA file
    Output:     blastn TXT output file    
    """

    stdout = os.path.join(blast_path,f"{input_args.array}.out")
    stderr = os.path.join(blast_path,f"{input_args.array}.err")

    cmd1 = f"singularity exec -B /nfs:/nfs {input_args.singularity} blastn -db {input_args.its_db} -query {filenames.its_fasta} -num_threads {input_args.cpus} -outfmt \"6 qacc sscinames sacc pident bitscore evalue stitle\" -max_target_seqs 5 -out {filenames.blast_out}"
    try:
        lib.print_n("BLASTn of ITS sequence against ITS DB")
        lib.execute_shell(cmd1,stdout,stderr)
    except subprocess.CalledProcessError as e:
        print(e.returncode)
        print(e.output)

def parse_blastp(results_path,array):
    """
    Parses BLASTp output files and appends the extracted information to a pandas DataFrame.

    Input:      BLASTp output file from blastp() function [path]
    Output:     pandas DataFrame object
    """

    columns = ["Query","Species","Accession","Identity","Bitscore","E-value","Subject Name"]

    if os.path.exists(results_path):
        lib.print_n(f"Parsing BLASTp results for {results_path}")
        try:
            df = pd.read_csv(results_path, sep="\t", header=None)
            df.columns = columns
            df = df.head(1)
            df["ID"] = array
            return df
        except pd.errors.EmptyDataError:
            print(f"No BLASTp results obtained for {results_path}")
    else:
        print(f"The BLASTp for {results_path} was not completed successfully")
        
def quast(input_args,filenames,output_path):
    """
    Runs `quast` genome analysis software on an assembly.

    Input:      assembly FASTA file and GFF file (optional)
    Output:     QUAST output folder
    """

    stdout = os.path.join(output_path,f"{input_args.array}.out")
    stderr = os.path.join(output_path,f"{input_args.array}.err")

    if os.path.exists(filenames.funannotate_gff) == True and os.stat(filenames.funannotate_gff).st_size > 0:
        lib.print_n(f"Funannotate annotations GFF file exists, analysing {filenames.funannotate_gff} and {filenames.assembly_fasta} with quast")
        cmd1 = ["singularity","exec","-B","/nfs:/nfs",input_args.singularity,"quast.py",filenames.assembly_fasta,"-o",output_path,"-g",filenames.funannotate_gff,"-t",input_args.cpus,"--fungus","-L","-b"]
    else:
        lib.print_n(f"Funannotate annotations GFF file does not exist, analysing {filenames.assembly_fasta} with quast")
        cmd1 = ["singularity","exec","-B","/nfs:/nfs",input_args.singularity,"quast.py",filenames.assembly_fasta,"-o",output_path,"-t",input_args.cpus,"--fungus","-L","-b"]

    try:
        lib.execute(cmd1,stdout,stderr)
        lib.file_exists_exit(filenames.quast_report,"Quast successfully analysed the assembly!","Quast failed... check the logs and your inputs")
    except subprocess.CalledProcessError as e:
        print(e.returncode)
        print(e.output) 

def parse_quast(filenames,array):
    """
    Parses `quast` output reports, returning a pandas dataframe.

    Input:      quast TSV report
    Output:     pandas DataFrame object
    """

    lib.print_n(f"Parsing Quast results for {filenames.quast_report}")
    try:
        df = pd.read_csv(filenames.quast_report, sep="\t", header=0)
        df["ID"] = array
        return df
    except pd.errors.EmptyDataError:
        lib.print_e(f"The QUAST report {filenames.quast_report} was not parsed successfully")

def antismash(input_args,filenames,antismash_out):
    """
    Runs `antismash` BGC analysis on assembly file.
        - antismash will have a cry about files already being present in the output folder, 
          so STDOUT and STDERR files found in parent dir, rather than 'antismash_out'.
    
    Input:      annotated assembly GBK file
    Output:     antismash output folder
    """
    
    stdout = os.path.join(input_args.directory_new,f"{input_args.array}_antismash.out")
    stderr = os.path.join(input_args.directory_new,f"{input_args.array}_antismash.err")

    cmd1 = ["singularity","exec","-B","/nfs:/nfs",input_args.singularity_antismash,"antismash","--cpus",input_args.cpus,"--taxon","fungi","--fullhmmer","--genefinding-tool","glimmerhmm","--cc-mibig","--smcog-trees","--cb-general","--cb-subclusters","--cb-knownclusters","--minlength","5000","--output-dir",antismash_out,filenames.antismash_assembly]

    try:
        lib.print_n(f"Predicting gene clusters from {filenames.antismash_assembly} with antiSMASH")
        lib.execute(cmd1,stdout,stderr)
        lib.file_exists_exit(filenames.antismash_index,"antiSMASH successfully analysed the assembly!","antiSMASH failed... check the logs and your inputs")
    except subprocess.CalledProcessError as e:
        print(e.returncode)
        print(e.output)

def main(input_args,filenames):

    lib.print_h("Initializing \'post analysis\' module...")
    post_start_time = datetime.datetime.now()
    quast_text = "  _______\n___/ Quast \____________________________________________________________________"
    lib.print_t(quast_text)

    quast_path = os.path.join(input_args.directory_new,"quast")
    filenames.quast_report = os.path.join(quast_path,"transposed_report.tsv")
    lib.make_path(quast_path)
    if lib.file_exists(filenames.quast_report,"Quast already analysed the assembly!","") is False:
        quast(input_args,filenames,quast_path)
        quast_df = parse_quast(filenames.quast_report,input_args.array)

    if len(input_args.its_db) > 0:
            
        tax_lookup_text = "  ___________________\n___/  Taxonomy Lookup  \________________________________________________________"        
        lib.print_t(tax_lookup_text)

        itsx_path = os.path.join(input_args.directory_new,"ITSx")
        lib.make_path(itsx_path)
        blastn_path = os.path.join(input_args.directory_new,"blastn")
        filenames.blast_out = os.path.join(blastn_path,f"{input_args.array}_its_blastn.out")
        its_full = os.path.join(itsx_path,f"{input_args.array}.full.fasta")
        its_2 = os.path.join(itsx_path,f"{input_args.array}.ITS2.fasta")
        its_1 = os.path.join(itsx_path,f"{input_args.array}.ITS1.fasta")
        if lib.file_exists(its_full,"ITSx already extracted the ITS sequence! Skipping","") is False:
            #itsx(input_args,filenames,itsx_path)
            lib.make_path(blastn_path)
            if lib.file_exists(its_full,"ITSx successfully extracted the full ITS sequence!","ITSx failed to extract the full ITS sequence. Checking for ITS2") == True:
                filenames.its_fasta = its_full
                blastn(input_args,filenames,blastn_path)
                if lib.file_exists(filenames.blast_out,f"BLASTn successfully searched {filenames.its_fasta} against the ITS_RefSeq database!",f"BLASTn failed to search {filenames.its_fasta}... check the logs and your inputs") is True:
                    its_df = parse_blastp(filenames.blast_out,input_args.array)
            elif lib.file_exists(its_2,"ITSx successfully extracted the ITS2 sequence!","ITSx failed to extract the ITS2 sequence. Checking for ITS1") is True:
                filenames.its_fasta = its_2
                blastn(input_args,filenames,blastn_path)
                if lib.file_exists(filenames.blast_out,f"BLASTn successfully searched {filenames.its_fasta} against the ITS_RefSeq database!",f"BLASTn failed to search {filenames.its_fasta}... check the logs and your inputs") is True:
                    its_df = parse_blastp(filenames.blast_out,input_args.array)
            elif lib.file_exists(its_1,"ITSx successfully extracted the ITS1 sequence!","ITSx failed... check the logs and your inputs") is True:
                filenames.its_fasta = its_1
                blastn(input_args,filenames,blastn_path)
                if lib.file_exists(filenames.blast_out,f"BLASTn successfully searched {filenames.its_fasta} against the ITS_RefSeq database!",f"BLASTn failed to search {filenames.its_fasta}... check the logs and your inputs") is True:
                    its_df = parse_blastp(filenames.blast_out,input_args.array)
    else:
        lib.print_h("Skipping taxonomy lookup of isolate")

    if input_args.antismash > 0:

        bgc_text = "  ____________\n___/ BGC search \_______________________________________________________________"
        lib.print_t(bgc_text)

        antismash_out = os.path.join(input_args.directory_new,"antismash")

        filenames.antismash_index = os.path.join(antismash_out,"index.html")
        if lib.file_exists(filenames.funannotate_gbk,f"Using {filenames.funannotate_gbk} as input for antiSMASH", \
            f"Using {filenames.assembly_fasta} as input for antiSMASH") == True:
            filenames.antismash_assembly = filenames.funannotate_gbk
        else:
            filenames.antismash_assembly = filenames.assembly_fasta
        if lib.file_exists(filenames.antismash_index,"antiSMASH already analysed the assembly!","") is False:
            if os.path.exists(antismash_out):
                lib.print_n("Removing existing antiSMASH output directory")
                shutil.rmtree(antismash_out)
            lib.make_path(antismash_out)
            antismash(input_args,filenames,antismash_out)
    else:
        lib.print_n("Skipping BGC search with anitSMASH. Use '--antismash' as a script argument if you would like to do this.")

    lib.print_h("Generating final dataframe")
    final_df = pd.DataFrame()
    if len(quast_df) > 0:
        final_df = quast_df
    if len(its_df) > 0:
        final_df = final_df.merge(its_df, on="ID")

    with pd.option_context('display.max_rows', None, 'display.max_columns', None):
        lib.print_n(final_df)
    
    results_path = os.path.join(input_args.directory_new,"results")
    lib.make_path(results_path)
    filenames.csv_output = os.path.join(results_path,f"{input_args.array}_results.csv")
    final_df.to_csv(filenames.csv_output, index=False)
    lib.file_exists_exit(filenames.csv_output,"Final results report successfully generated","Final results report was not generated")
    lib.print_h(f"post analysis module completed in {datetime.datetime.now() - post_start_time}")

if __name__ == '__main__':
    main()