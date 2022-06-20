import os
import sys, argparse
import library as lib

"""
This script will install the databases required to operate the Fungiflow pipeline.

The estimated time to download and prepare all the databases will take 2 days with 8 cpus

If you do not requre all the databases (e.g., no NCBI-nt because you don't want to use the blobplot module),
you can install a custom set of databases by adding a set of database strings to the `--databases` parameter as below:
    -all        Installs all databases
    -kraken2    Installs the Kraken2 database for the main module (optional)
    -ncbi-nt    Installs the NCBI-nt database for the blobplot module (optional)
    -ncbi-its   Installs the NCBI-ITS-refseq database for the post_analysis module (optional)

    To install only the kraken2 and ncbi-its databases, for example, use the following:
        `--databases "kraken2,ncbi-its"`
"""
def get_args():
    """Parse command line arguments"""
    
    try:
        parser = argparse.ArgumentParser(
            description="Fungiflow - the automated eukaryotic genomic pipeline for fungi.")
        parser.add_argument('-d', '--directory', action='store',
                            help='Database directory path.', type=str, default="databases")  
        parser.add_argument('-db', '--databases', action='store',
                            help='Databases to install. Accepted arguments are "all" or any combo of "kraken2", \
                            "ncbi-its","ncbi-nt".', type=str, required=True)   
        parser.add_argument('-c', '--cpus', action='store',
                            help='Number of threads to use', type=str, required=True)
        parser.add_argument('-m', '--mem', action='store',
                            help='Amount of memory to use (in GB)', type=str, required=True)
     
    except argparse.ArgumentError:
        lib.print_e("An exception occurred with argument parsing. Check your inputs.")
        exit(1)

    def __str__(self):
        return  str(self.__class__) + '\n' + '\n'.join((str(item) + ' = ' + str(self.__dict__[item]) for item in sorted(self.__dict__)))

    return parser.parse_args()

def install_kraken2_db(input_args,database_path):
    """
    Downloads the Kraken2 standard database.
    """

    stdout = "kraken2_db.out"
    stderr = "kraken2_db.err"

    cmd1 = ["singularity","exec","-B","/nfs:/nfs",input_args.singularity,"kraken2-build","--use-ftp","-i","--standard","--threads",input_args.cpus,"--db",database_path]
    cmd2 = ["singularity","exec","-B","/nfs:/nfs",input_args.singularity,"kraken2-build","--build","--threads",input_args.cpus,"--db",database_path]
    cmd3 = ["singularity","exec","-B","/nfs:/nfs",input_args.singularity,"kraken2-build","--clean","--threads",input_args.cpus,"--db",database_path]
    cmd4 = ["singularity","exec","-B","/nfs:/nfs",input_args.singularity,"kraken2-build","----download-library","human","--no-masking","--threads",input_args.cpus,"--db",database_path]
    cmd5 = ["singularity","exec","-B","/nfs:/nfs",input_args.singularity,"kraken2-build","----download-library","UniVec","--no-masking","--threads",input_args.cpus,"--db",database_path]

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
    
    cmd1 = ["singularity","exec","-B","/nfs:/nfs",input_args.singularity,"update_blastdb.pl","--passive","--decompress","ITS_RefSeq_Fungi"]
    cmd2 = ["singularity","exec","-B","/nfs:/nfs",input_args.singularity,"update_blastdb.pl","taxdb"]
    cmd3 = ["tar","-xzf","taxdb.tar.gz"]

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
    
    cmd1 = ["singularity","exec","-B","/nfs:/nfs",input_args.singularity,"update_blastdb.pl","--passive","--decompress","nt" ]
    cmd2 = ["singularity","exec","-B","/nfs:/nfs",input_args.singularity,"update_blastdb.pl","taxdb"]
    cmd3 = ["tar","-xzf","taxdb.tar.gz"]

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

def inspect(databases):
    """
    Confirms that the databases have been installed and are correct.
    """

    stdout = "inspect_dbs.out"
    stderr = "inspect_dbs.err"

    cmds = []
    if "kraken2" in databases:
        kraken2 = ["kraken2-inspect",]
        cmds.append(kraken2)
    if "its" in databases:
        its = ["its",]
        cmds.append(its)
    if "nt" in databases:
        nt = ["blastdbcmd",]
        cmds.append(nt)

    try:
        lib.execute(cmds,stdout,stderr)
    except subprocess.CalledProcessError as e:
        print(e.returncode)
        print(e.output)

def main(input_args):

    if input_args is None:
        input_args = get_args()

    args = get_args()
    lib.print_tu("\n⁂⁂⁂⁂⁂⁂⁂⁂ Install Script Begins ⁂⁂⁂⁂⁂⁂⁂⁂\n")
    start_time = datetime.datetime.now()
    os.chdir(args.directory)
    lib.print_n(args)
    databases_path = os.path.join("databases")



    # need to check this code block to see if what dbs are printed with a given input
    dbs = ["kraken2","ncbi-its","ncbi-nt"]
    if args.databases is not "all":
        for i in args.databases:
            if i not in dbs:
                lib.print_e(f"{i} is not a valid database option. Valid options are {dbs}")
                exit()
        dbs = args.databases
    
    if "kraken2" in dbs:
        kraken2_time = datetime.datetime.now()
        print("Downloading and installing Kraken2 database...")
        kraken2_path = os.path.join(databases_path,"kraken2")
        install_kraken2_db(args,kraken2_path)
        lib.print_h(f"Kraken2 database installed in {datetime.datetime.now() - kraken2_time}")
    if "ncbi-its" in dbs:
        its_time = datetime.datetime.now()
        print("Downloading and installing NCBI-ITSrefseq database...")
        its_path = os.path.join(databases_path,"ncbi-its")
        os.chdir(its_path)
        install_ncbi_its(args)
        lib.print_h(f"NCBI-ITSrefseq database installed in {datetime.datetime.now() - its_time}")
        os.chdir("..")
    if "ncbi-nt" in dbs:
        nt_time = datetime.datetime.now()
        print("Downloading and installing NCBI-nt database...")
        nt_path = os.path.join(databases_path,"ncbi-nt")
        os.chdir(nt_path)
        install_ncbi_nt(args)
        lib.print_h(f"NCBI-nt database installed in {datetime.datetime.now() - nt_time}")
        os.chdir("..")
    
    inspect(dbs,paths)

    lib.print_h(f"All databases installed in {datetime.datetime.now() - start_time}")
    lib.print_tu("\n⁂⁂⁂⁂⁂⁂⁂⁂ Script Finished ⁂⁂⁂⁂⁂⁂⁂⁂\n")

if __name__ == '__main__':
    main()