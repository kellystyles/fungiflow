#!/usr/bin/env python3
import os
import argparse
import subprocess
import datetime
import library as lib

############## Need to add eggnog and interproscan dbs

"""
This script will install the databases required to operate the Fungiflow pipeline.

The estimated time to download and prepare all the databases will take 2 days with 8 cpus

If you do not requre all the databases (e.g., no NCBI-nt because you don't want to use the blobplot module),
you can install a custom set of databases by adding a set of database strings to the `--databases` parameter as below:
    -all        Installs all databases
    -kraken2    Installs the Kraken2 database for the main module (optional)
    -ncbi-nt    Installs the NCBI-nt database for the blobplot module (optional)
    -ncbi-its   Installs the NCBI-ITS-refseq database for the post_analysis module (optional)
    -eggnog     Installs the EggNOG database for functional annotation (optional)

    To install only the kraken2 and ncbi-its databases, for example, use the following:
        `--databases "kraken2,ncbi-its"`
"""
def get_args():
    """Parse command line arguments"""
    
    try:
        parser = argparse.ArgumentParser(
            description="Fungiflow - the automated eukaryotic genomic pipeline for fungi.")
        parser.add_argument('-d', '--directory', action='store',
                            help='Database directory path.', type=str, default="databases", required=False)  
        parser.add_argument('-db', '--databases', action='store',
                            help='Databases to install. Accepted arguments are "all" or any combo of "kraken2", \
                            "ncbi-its","ncbi-nt", or "eggnog.', type=str, default="all", required=False)   
        parser.add_argument('-c', '--cpus', action='store',
                            help='Number of threads to use', type=str, default="4", required=False)
        parser.add_argument('-m', '--mem', action='store',
                            help='Amount of memory to use (in GB)', type=str, default="8", required=False)
        parser.add_argument('-s', '--singularity_fungiflow', action='store',
                            help='Path to Fungiflow Singularity container', type=str, required=False)                            
        parser.add_argument('-sfun', '--singularity_funannotate', action='store',
                            help='Path to Funannotate Singularity container', type=str, required=False)        
    except argparse.ArgumentError:
        lib.print_e("An exception occurred with argument parsing. Check your inputs.")
        exit(1)

    def __str__(self):
        return  str(self.__class__) + '\n' + '\n'.join((str(item) + ' = ' + str(self.__dict__[item]) for item in sorted(self.__dict__)))

    return parser.parse_args()

def install_kraken2_db(input_args,database_path):
    """
    Downloads the Kraken2 standard database.
    Will first attempt to download the standard preformatted database (cmd1), else will attempt to download this same database library by library (cmds 4 & 5).
    The library will then be built and cleaned (cmds 2 & 3)
    """

    stdout = "kraken2_db.out"
    stderr = "kraken2_db.err"

    cmd1 = ["kraken2-build","--use-ftp","--standard","--threads",input_args.cpus,"--db",database_path]
    cmd2 = ["kraken2-build","--build","--threads",input_args.cpus,"--db",database_path]
    cmd3 = ["kraken2-build","--clean","--threads",input_args.cpus,"--db",database_path]
    cmd4 = ["kraken2-build","----download-library","human","--no-masking","--threads",input_args.cpus,"--db",database_path]
    cmd5 = ["kraken2-build","----download-library","UniVec","--no-masking","--threads",input_args.cpus,"--db",database_path]
    if len(input_args.singularity_fungiflow) > 0: cmd1 = input_args.singularity1 + cmd1
    if len(input_args.singularity_fungiflow) > 0: cmd2 = input_args.singularity1 + cmd2
    if len(input_args.singularity_fungiflow) > 0: cmd3 = input_args.singularity1 + cmd3
    if len(input_args.singularity_fungiflow) > 0: cmd4 = input_args.singularity1 + cmd4
    if len(input_args.singularity_fungiflow) > 0: cmd5 = input_args.singularity1 + cmd5

    print(" ".join(cmd1))
    print(" ".join(cmd2))
    print(" ".join(cmd3))
    print(" ".join(cmd4))
    print(" ".join(cmd5))

    try:
        lib.execute(cmd1,stdout,stderr)
    except subprocess.CalledProcessError as e:
        print(e.returncode)
        print(e.output)
        lib.execute(cmd4,stdout,stderr)
        lib.execute(cmd5,stdout,stderr)
    try:
        lib.execute(cmd2,stdout,stderr)
        lib.execute(cmd3,stdout,stderr)
    except subprocess.CalledProcessError as e:
        print(e.returncode)
        print(e.output)

def install_ncbi_its(input_args):
    """
    Downloads the NCBI ITS fungi refseq preformatted BLASTn database.
    """

    stdout = "ncbi-its_db.out"
    stderr = "ncbi-its_db.err"
    
    cmd1 = ["update_blastdb.pl","--passive","--decompress","ITS_RefSeq_Fungi"]
    cmd2 = ["update_blastdb.pl","taxdb"]
    cmd3 = ["tar","-xzf","taxdb.tar.gz"]
    if len(input_args.singularity_fungiflow) > 0: cmd1 = input_args.singularity1 + cmd1
    if len(input_args.singularity_fungiflow) > 0: cmd2 = input_args.singularity1 + cmd2
    print(" ".join(cmd1))
    print(" ".join(cmd2))
    print(" ".join(cmd3))

    try:
        lib.execute(cmd1,stdout,stderr)
        lib.execute(cmd2,stdout,stderr)
        lib.execute(cmd3,stdout,stderr)
    except subprocess.CalledProcessError as e:
        print(e.returncode)
        print(e.output)

def install_ncbi_nt(input_args):
    """
    Downloads the NCBI nt preformatted BLASTn database.
    """

    stdout = "ncbi-nt_db.out"
    stderr = "ncbi-nt_db.err"
    
    cmd1 = ["update_blastdb.pl","--passive","--decompress","nt" ]
    cmd2 = ["update_blastdb.pl","taxdb"]
    cmd3 = ["tar","-xzf","taxdb.tar.gz"]
    if len(input_args.singularity_fungiflow) > 0: cmd1 = input_args.singularity1 + cmd1
    if len(input_args.singularity_fungiflow) > 0: cmd2 = input_args.singularity1 + cmd2

    print(" ".join(cmd1))
    print(" ".join(cmd2))
    print(" ".join(cmd3))

    try:
        lib.execute(cmd1,stdout,stderr)
        lib.execute(cmd2,stdout,stderr)
        lib.execute(cmd3,stdout,stderr)
    except subprocess.CalledProcessError as e:
        print(e.returncode)
        print(e.output)

def install_eggnog(input_args, database_path):
    """
    Downloads the EggNOG database.
    """

    stdout = "eggnog_db.out"
    stderr = "eggnog_db.err"
    
    cmd1 = ["export", f"EGGNOG_DATA_DIR={database_path}"]
    cmd2 = ["create_dbs.py", "-m", "diamond", "--dbname", "fungi", "--taxa", "Fungi"]
    if len(input_args.singularity_funannotate) > 0: cmd2 = input_args.singularity2 + cmd2
    print(" ".join(cmd1))
    print(" ".join(cmd2))

    try:
        lib.execute(cmd1,stdout,stderr)
        lib.execute(cmd2,stdout,stderr)
    except subprocess.CalledProcessError as e:
        print(e.returncode)
        print(e.output)

def main():

    args = get_args()
    lib.print_tu("\n⁂⁂⁂⁂⁂⁂⁂⁂ Install Script Begins ⁂⁂⁂⁂⁂⁂⁂⁂\n")
    start_time = datetime.datetime.now()
    os.chdir(args.directory)
    lib.print_n(args)
    databases_path = os.path.join(args.directory, "databases")

    if len(args.singularity_fungiflow) > 0:
        args.singularity1 = ["singularity", "exec", args.singularity_fungiflow]
    if len(args.singularity_funannotate) > 0:
        args.singularity2 = ["singularity", "exec", args.singularity_funannotate]

    # need to check this code block to see if what dbs are printed with a given input
    dbs = ["kraken2","ncbi-its","ncbi-nt","eggnog"]
    if args.databases != "all":
        for i in args.databases:
            if i not in dbs:
                lib.print_e(f"{i} is not a valid database option. Valid options are {dbs}")
                exit()
        dbs = args.databases
    
    if "kraken2" in dbs:
        kraken2_time = datetime.datetime.now()
        print("Downloading and installing Kraken2 database...")
        kraken2_path = os.path.join(databases_path,"kraken2")
        os.makedirs(kraken2_path)     
        install_kraken2_db(args,kraken2_path)
        lib.print_h(f"Kraken2 database installed in {datetime.datetime.now() - kraken2_time}")
    if "ncbi-its" in dbs:
        its_time = datetime.datetime.now()
        print("Downloading and installing NCBI-ITSrefseq database...")
        its_path = os.path.join(databases_path,"ncbi-its")
        os.makedirs(its_path)        
        os.chdir(its_path)
        install_ncbi_its(args)
        lib.print_h(f"NCBI-ITSrefseq database installed in {datetime.datetime.now() - its_time}")
        os.chdir(databases_path)
    if "ncbi-nt" in dbs:
        nt_time = datetime.datetime.now()
        print("Downloading and installing NCBI-nt database...")
        nt_path = os.path.join(databases_path,"ncbi-nt")
        os.makedirs(nt_path)
        os.chdir(nt_path)
        install_ncbi_nt(args)
        lib.print_h(f"NCBI-nt database installed in {datetime.datetime.now() - nt_time}")
        os.chdir(databases_path)
    if "eggnog" in dbs:
        eggnog_time = datetime.datetime.now()
        print("Downloading and installing EggNOG database...")
        eggnog_path = os.path.join(databases_path,"eggnog")
        os.makedirs(eggnog_path)
        os.chdir(eggnog_path)
        install_eggnog(args, eggnog_path)
        lib.print_h(f"NCBI-nt database intalled in {datetime.datetime.now() - eggnog_time}")
        os.chdir(databases_path)
    lib.check_databases(args)

    lib.print_h(f"All databases installed in {datetime.datetime.now() - start_time}")
    lib.print_tu("\n⁂⁂⁂⁂⁂⁂⁂⁂ Script Finished ⁂⁂⁂⁂⁂⁂⁂⁂\n")

if __name__ == '__main__':
    main()