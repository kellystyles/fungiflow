#!/usr/bin/env python3
import os
import argparse
import subprocess
import datetime
import requests
import tarfile
import library as lib

"""
This script will install the databases required to operate the Fungiflow pipeline.

The estimated time to download and prepare all the databases will take 2 days with 8 cpus

If you do not requre all the databases (e.g., no NCBI-nt because you don't want to use the blobplot module),
you can install a custom set of databases by adding a set of database strings to the `--databases` parameter as below:
    -all        Installs all databases
    -kraken2    Installs the standard Kraken2 database for the main module (optional)
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

def download_file(url, save_path):
    response = requests.get(url, stream=True)
    response.raise_for_status()
    
    with open(save_path, 'wb') as file:
        for chunk in response.iter_content(chunk_size=8192):
            file.write(chunk)

def extract_tar_gz(file_path, extract_path):
    with tarfile.open(file_path, 'r:gz') as tar:
        tar.extractall(extract_path)

def install_kraken2_db(database_path):
    """
    Downloads the Kraken2 standard database.
    Will first attempt to download the standard preformatted database (cmd1), else will attempt to download this same database library by library (cmds 4 & 5).
    The library will then be built and cleaned (cmds 2 & 3)
    """

    # GitHub URL
    url = 'https://github.com/BenLangmead/aws-indexes/blob/master/docs/k2.md'

    # Send a GET request to fetch the HTML content
    response = requests.get(url)
    html_content = response.text

    # parse data and get links
    parsed = html_content.split("\\")
    stripped_list = [item.strip('"') for item in parsed]
    links = []
    for i in stripped_list:
        if "k2_standard_16" in i:
            links.append(i)
    latest = links[0]
    latest_output = latest.split("/")[-1]

    try:
        download_file(latest, latest_output)
        extract_tar_gz(latest_output, database_path)
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
    if len(input_args.singularity_fungiflow) > 0: cmd1 = input_args.singularity1 + cmd1
    if len(input_args.singularity_fungiflow) > 0: cmd2 = input_args.singularity1 + cmd2

    try:
        lib.execute(cmd1,stdout,stderr)
        lib.execute(cmd2,stdout,stderr)
        extract_tar_gz("taxdb.tar.gz", ".")
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

    cmd1 = ["download_eggnog_data.py", "-y", "--data_dir", "."]    
    cmd2 = ["create_dbs.py", "-m", "diamond", "--dbname", "fungi", "--taxa", "Fungi", "-y", "--data_dir", "."]
    if len(input_args.singularity_funannotate) > 0: cmd1 = input_args.singularity2 + cmd1
    if len(input_args.singularity_funannotate) > 0: cmd2 = input_args.singularity2 + cmd2

    try:
        lib.execute(cmd1,stdout,stderr)
        lib.execute(cmd2,stdout,stderr)
        print(f"\nExport the GENEMARK_PATH to your environment like so:\n'export GENEMARK_PATH=${database_path}'\n")
    except subprocess.CalledProcessError as e:
        print(e.returncode)
        print(e.output)

def check_databases(input_databases, databases_path):
    """
    Will check key files in each database in 'input_databases' at 
    'databases_path' and will return add the database name to either a 'passed'
    or 'failed' list.
    """
    
    failed = []
    passed = []
    #print(input_databases)

    # define database paths
    kraken2_path = os.path.join(databases_path,"kraken2")
    its_path = os.path.join(databases_path,"ncbi-its")
    nt_path = os.path.join(databases_path,"ncbi-nt")    
    eggnog_path = os.path.join(databases_path,"eggnog")

    try:
        if "ncbi-nt" in input_databases:
            ncbi_nt = os.path.join(nt_path, "nt.00.nhd")
            taxdb = os.path.join(nt_path, "taxdb.bti")
            if lib.file_exists_list([ncbi_nt, taxdb], "", "") is False:
                failed.append("ncbi-nt")    
            elif lib.file_exists_list([ncbi_nt, taxdb], "", "") is True:
                passed.append("ncbi-nt")        
                #print(nt_path)
    except UnboundLocalError:
        pass
    try:
        if "ncbi-its" in input_databases:
            ncbi_its = os.path.join(its_path, "ITS_RefSeq_Fungi.nsq")
            taxdb = os.path.join(its_path, "taxdb.bti")
            if lib.file_exists_list([ncbi_its, taxdb], "", "") is False:
                failed.append("ncbi-its")
            elif lib.file_exists_list([ncbi_its, taxdb], "", "") is True:
                passed.append("ncbi-its") 
                #print(its_path)
    except UnboundLocalError:
        pass
    try:
        if "kraken2" in input_databases:
            hash = os.path.join(kraken2_path, "hash.k2d")
            opts = os.path.join(kraken2_path, "opts.k2d")
            taxo = os.path.join(kraken2_path, "taxo.k2d")
            if lib.file_exists_list([hash, opts, taxo], "", "") is False:
                failed.append("kraken2")
            elif lib.file_exists_list([hash, opts, taxo], "", "") is True:
                passed.append("kraken2")  
                #print(kraken2_path)
    except UnboundLocalError:
        pass
    try:
        if "eggnog" in input_databases:
            dmnd = os.path.join(eggnog_path, "fungi.dmnd")
            taxa = os.path.join(eggnog_path, "eggnog.taxa.db")
            main = os.path.join(eggnog_path, "eggnog.db")
            prots = os.path.join(eggnog_path, "e5.proteomes.faa")
            fungi_dmnd = os.path.join(eggnog_path, "fungi.dmnd")
            if lib.file_exists_list([dmnd, taxa, main, prots, fungi_dmnd], "", "") is False:
                failed.append("eggnog")
            elif lib.file_exists_list([dmnd, taxa, main, prots], "", "") is True:
                passed.append("eggnog") 
                #print(eggnog_path)
    except UnboundLocalError:
        pass

    return passed, failed

def main():

    # get arguments
    args = get_args()
    lib.print_t("\n⁂⁂⁂⁂⁂⁂⁂⁂ Install Script Begins ⁂⁂⁂⁂⁂⁂⁂⁂\n")
    start_time = datetime.datetime.now()

    # define the singularity container paths
    if args.singularity_fungiflow is not None:
        args.singularity1 = ["singularity", "exec", os.path.abspath(os.path.join(args.singularity_fungiflow))]
    if args.singularity_funannotate is not None:
        args.singularity2 = ["singularity", "exec", os.path.abspath(os.path.join(args.singularity_funannotate))]

    # create/move to the databases directory
    if not os.path.exists(args.directory):
        os.makedirs(args.directory)
    #print(args)
    databases_path = os.path.abspath(os.path.join(args.directory))
    print(f"Database path: {databases_path}")
    os.chdir(databases_path)

    # need to check this code block to see if what dbs are printed with a given input
    dbs = ["kraken2","ncbi-its","ncbi-nt","eggnog"]
    if args.databases != "all":
        input_databases = args.databases.split(',')
        for i in input_databases:
            if i.strip() not in dbs:
                lib.print_e(f"{i.strip()} is not a valid database option. Valid options are {dbs}")
                exit()
        dbs = args.databases
    else:
        input_databases = args.databases

    # define database paths
    kraken2_path = os.path.join(databases_path,"kraken2")
    its_path = os.path.join(databases_path,"ncbi-its")
    nt_path = os.path.join(databases_path,"ncbi-nt")    
    eggnog_path = os.path.join(databases_path,"eggnog")

    # check what databases exist
    passed, failed = check_databases(input_databases, databases_path)
    #print(passed, failed)
    
    # install databases
    if len(passed) > 0:
        print("\nThe following databases exist:")
        for db in passed:
            print(f" - {db}, skipping installation...")
    if len(failed) > 0:
        print("\nInstalling the following databases:")
        for db in failed:
            print(f" - {db}")
    print("\n")    
    if "kraken2" in dbs and "kraken2" in failed:
        kraken2_time = datetime.datetime.now()
        print("Downloading and installing kraken2 database...")
        if not os.path.exists(kraken2_path):
            os.makedirs(kraken2_path)    
        install_kraken2_db(kraken2_path)
        lib.print_h(f"kraken2 database installed in {datetime.datetime.now() - kraken2_time}")
    if "ncbi-its" in dbs and "ncbi-its" in failed:
        its_time = datetime.datetime.now()
        print("Downloading and installing ncbi-its database...")
        if not os.path.exists(its_path):
            os.makedirs(its_path)       
        os.chdir(its_path)
        install_ncbi_its(args)
        lib.print_h(f"ncbi-its database installed in {datetime.datetime.now() - its_time}")
        os.chdir(databases_path)
    if "ncbi-nt" in dbs and "ncbi-nt" in failed:
        nt_time = datetime.datetime.now()
        print("Downloading and installing ncbi-nt database...")
        if not os.path.exists(nt_path):
            os.makedirs(nt_path) 
        os.chdir(nt_path)
        install_ncbi_nt(args)
        lib.print_h(f"ncbi-nt database installed in {datetime.datetime.now() - nt_time}")
        os.chdir(databases_path)
    if "eggnog" in dbs and "eggnog" in failed:
        eggnog_time = datetime.datetime.now()
        print("Downloading and installing the eggnog database...")
        if not os.path.exists(eggnog_path):
            os.makedirs(eggnog_path) 
        os.chdir(eggnog_path)
        install_eggnog(args, eggnog_path)
        lib.print_h(f"eggnog database intalled in {datetime.datetime.now() - eggnog_time}")
        os.chdir(databases_path)
    
    # check all databases were installed as expected
    passed, failed = check_databases(input_databases, databases_path)
    
    # states which database(s) were installed correctly.
    if len(passed) > 0:
        print("The following databases were successfully installed/existed:")
        for db in passed:
            print(f" - {db} is installed!")
    if len(failed) > 0:
        lib.print_e("\nThere are issues with the following databases:")
        for i in failed:
            print(f" - {i}")
        lib.print_e("Please check the installation path.")
        lib.print_e("Try removing the suspect database directories and run `install.py` again.")
    
    if len(failed) == 0:
        lib.print_h(f"\nAll databases installed in {datetime.datetime.now() - start_time}")
    lib.print_t("\n⁂⁂⁂⁂⁂⁂⁂⁂ Script Finished ⁂⁂⁂⁂⁂⁂⁂⁂\n")

if __name__ == '__main__':
    main()