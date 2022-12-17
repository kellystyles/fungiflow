import os
import sys
import subprocess
import datetime

"""
Helper functions for the Fungiflow pipeline.
"""

# Creating the Class for the workflow filepaths
class Files:
    """
    This Class object will take as input the short read file paths and from there can be used to store all
    other pipeline filenames.
    It also will initialize lists to append Singularity commands to the front of relevant BASH commands if
    Singularity container paths are supplied.
    This Class has two printer functions, one for printing the current filenames stored in the Class object,
    and one for printing all the files from this Class that exist and have a size >0.
    """
    def __init__(self, input_args):
        # defines paths for short reads
        self.shortf = os.path.abspath(os.path.join(input_args.directory, input_args.illumina_f))
        self.shortr = os.path.abspath(os.path.join(input_args.directory, input_args.illumina_r))     
        # defines whether singularity commands are needed to execute the pipeline
        if input_args.singularity_image is not None: 
            self.singularity = ["singularity", "exec", input_args.singularity_image]
        else:
            self.singularity = ""
        if input_args.singularity_funannotate is not None: 
            self.funannotate = ["singularity", "exec", input_args.singularity_funannotate]
        else:
            self.funannotate = ""
        if input_args.singularity_antismash is not None: 
            self.antismash = ["singularity", "exec", input_args.singularity_antismash]
        else:
            self.antismash = ""
    
    # function to print all filenames currently in this Class object
    def __str__(self):
        return  str(self.__class__) + '\n' + '\n'.join((str(item) + ' = ' + str(self.__dict__[item]) for item in sorted(self.__dict__)))
    
    # printer function that will print each file in filenames that exists and has a size >0
    def printer(self):
        for k,  v in vars(self).items():
            try:
                if os.path.isfile(v) is True and os.stat(v).st_size > 0:
                    print(k, "\t", v)
            except TypeError:
                print(k," ".join(v))
                
def make_path(path):
    """
    Makes a path, regardless of if the parent dir exists or not.
    Functionally equivalent to `mkdir -p` in BASH
    """
    try:
        os.makedirs(path)
    except FileExistsError:
        pass

def file_exists(file, success_message, error_message):
    """
    Check if file exists and is not empty.
    Returns success message if the above is True.
    Returns error message if the above is False.
    """    
    
    if os.path.isfile(file) is True and os.stat(file).st_size > 0:
        if len(success_message) > 1:
            print_n(success_message)
        return True
    else:
        if len(error_message) > 1:
            print_e(error_message)
        return False

def file_exists_list(file_list, success_message, error_message):
    """
    Will check if a list of files exist.
    """
    # append boolean result of `file_exists` func to a dict for each file in list
    d = {}
    for i in file_list:
        bool = file_exists(i,"","")
        d[i] = bool
    # if 'False' in the dict, print error message, else print success message,
    # returning the appropriate bool result
    if False in d:
        if len(error_message) > 1:
            print_e(error_message)
        return False
    else:        
        if len(success_message) > 1:
            print_n(success_message)
        return True

def file_exists_exit(file, success_message, error_message):
    """
    Check if file exists and is not empty.
    Returns success message if the above is True.
    Returns error message if the above is False and will exit script.
    """
    if os.path.isfile(file) is True and os.stat(file).st_size > 0:
        print_n(success_message)
        return True
    else:
        print_e(error_message)
        sys.exit()

def file_exists_bool(filea, fileb, modifier, success_message, error_message):
    """
    Check if fileA exists and is not smaller than fileB by a modifier value.
    Returns error message if the above is not True.
    """
    try:
        if os.stat(filea).st_size > (os.stat(fileb).st_size * modifier):
            print_n(success_message)
            return True
        elif len(error_message) > 1:
            print_e(error_message)
    except FileNotFoundError:
        if len(error_message) > 1:
            print_e(error_message)
        return False

def execute(command, stdout, stderr):
    """
    Executes a Bash command and directs the stdout and stderr to separate files,  
    respectively.

    Parameters
    ----------
    command : comma separated list of strings (e.g.,  ["blastn", "-query", infile, "-out", outfile]).
    stdout : file to record standard out.
    stderr : file to record standard error.
    """
    try:
        print(" ".join(command))
        with open(stdout,  "wt") as out,  open(stderr,  "wt") as err:
            subprocess.run(command, stdout=out, stderr=err)
    except subprocess.CalledProcessError as e:
        print(e.returncode)
        print(e.output)   
        print(sys.exc_info()[1])

def execute_shell(command_list, stdout, stderr):
    """"
    Executes a list of bash commands, piping each command in the list after the first etc.
    i.e., command[1] will take the stdout of command[0] as stdin
    
    Commands must be in string format, with spaces included
    """
    
    # opens the stdout and stderr files for writing
    try:
        with open(stdout,  "wt") as out,  open(stderr,  "wt") as err:
            n = 2
            # If just a single command in the list will execute without PIPE
            if len(command_list) == 1:
                print(command_list[0])
                subprocess.run(command_list[0], stdout=out, stderr=err, shell=True)
            else:
                # takes first command in command_list as first command, notice no `input=` as with later processes
                # uses `exec` rather than subprocess.run to enable output to be piped to stdout
                first = f"proc1 = subprocess.run([\"{command_list[0]}\"], stdout=subprocess.PIPE, stderr=err, shell=True)"
                print(command_list[0])
                exec(first)
                # loops through remaining commands, excluding the last one
                if len(command_list[1:-1]) > 2:
                    for i in command_list[1:-1]:
                        # takes output from previous process as input, outputting to stdout via PIPE (eqv. to `|` in shell)
                        print(i)
                        nom = f"proc{n} = subprocess.run([\"{i}\"], input=proc{n-1}.stdout, stdout=subprocess.PIPE, stderr=err, shell=True)"
                        n += 1
                        exec(nom)
                # executing last command in command_list, this time writing to stdout, otherwise file won't be saved
                if len(command_list) > 1:
                    last = f"proc{n} = subprocess.run([\"{command_list[-1]}\"], input=proc{n-1}.stdout, stdout=out, stderr=err, shell=True)"
                    print(command_list[-1])
                    exec(last)
    except subprocess.CalledProcessError as e:
        print(e.returncode)
        print(e.output)   
        print(sys.exc_info()[1])     

def pickling(input_args, filenames):
    """
    Function that saves the class objects `input_args` and `filenames` to a TXT 
    file. 
    """
    inputs = open("input_args.txt",  'wb')
    pickle.dump(input_args, inputs)
    files = open("filenames.txt",  'wb')
    pickle.dump(filenames, files)
    inputs.close()
    files.close()

def unpickling(inputs_txt, files_txt):
    """
    Function that opens previously pickled class objects,  `input_args` and 
    `filenames`,  from TXT files.
    """
    inputs_f = open(inputs_txt,  'rb') 
    input_args = pickle.load(inputs_f)
    files_f = open(files_txt,  'rb') 
    filenames = pickle.load(files_f)

def check_databases(input_args):
    """
    Checks if the necessary databases are present based on input arguments.
    Then checks if the database exists and will assign the database to input_args.
    """
    # TRY statements to check what DB paths are supplied
    # In the case where the `--database_path` argument is not supploed, 
    # the supplied individual database paths will be used
    try: 
        if input_args.database_path is not None:
            if input_args.blobplot is True:
                blob_db = os.path.join(input_args.database_path, "NCBI_nt", "nt")
            if input_args.its is True:
                its_db = os.path.join(input_args.database_path, "fungi_ITS", "ITS_RefSeq_Fungi")
            if input_args.kraken2 is True:
                kraken2_db = os.path.join(input_args.database_path, "kraken2_std")
            if input_args.eggnog is True:
                eggnog_db = os.path.join(input_args.database_path, "eggnog")
        else:
            try:
                blob_db = input_args.blobplot_db
            except AttributeError:
                pass
            try:
                its_db = input_args.its_db
            except AttributeError:
                pass
            try: 
                kraken2_db = input_args.kraken2_db
            except AttributeError:
                pass
            try: 
                eggnog_db = input_args.eggnog_db
            except AttributeError:
                pass
    except IndexError:
        print("Please supply a value for `--database_path`")

    # Checks if key files are present in each database and assigns the main
    # `input_args` Class the individual databases.
    # Should the database files not exist, as determined by `file_exists_list`,
    # will add this database path to a list.
    c = []
    try:
        ncbi_nt = os.path.join(input_args.database_path, "NCBI_nt", "nt.00.nhd")
        taxdb = os.path.join(input_args.database_path, "NCBI_nt", "taxdb.bti")
        if file_exists_list([ncbi_nt, taxdb], "Blobplot DB is present", "Blobplot DB is not present. \
            Please supply a path via `--blobplot_db` or run install.py again.") is False:
            c.append(blob_db)            
        else:
            input_args.blobplot_db = blob_db
        print(blob_db)
    except UnboundLocalError:
        pass
    try:
        ncbi_its = os.path.join(input_args.database_path, "fungi_ITS", "ITS_RefSeq_Fungi.nsq")
        taxdb = os.path.join(input_args.database_path, "fungi_ITS", "taxdb.bti")
        if file_exists_list([ncbi_its, taxdb], "ITS DB is present", "ITS DB is not present. \
            Please supply a path via `--its_db` or run install.py again.") is False:
            c.append(its_db)
        else:
            input_args.its_db = its_db
        print(its_db)
    except UnboundLocalError:
        pass
    try:
        hash = os.path.join(kraken2_db, "hash.k2d")
        opts = os.path.join(kraken2_db, "opts.k2d")
        taxo = os.path.join(kraken2_db, "taxo.k2d")
        if file_exists_list([hash, opts, taxo], "kraken2 DB is present", "kraken2 DB is not present. \
            Please supply a path via `--kraken2_db` or run install.py again.") is False:
            c.append(kraken2_db)
        else:
            input_args.kraken2_db = kraken2_db    
        print(kraken2_db)
    except UnboundLocalError:
        pass
    try:
        dmnd = os.path.join(eggnog_db, "eggnog_proteins.dmnd")
        taxa = os.path.join(eggnog_db, "eggnog.taxa.db")
        main = os.path.join(eggnog_db, "eggnog.db")
        prots = os.path.join(eggnog_db, "e5.proteomes.faa")
        if file_exists_list([dmnd, taxa, main, prots], "eggnog DB is present", "eggnog DB is not present. \
            Please supply a path via `--eggnog_db` or run install.py again.") is False:
            c.append(eggnog_db)
        else:
            input_args.eggnog_db = eggnog_db    
        print(eggnog_db)
    except UnboundLocalError:
        pass
    # will exit program if one or more database(s) was missing files.
    if len(c) > 0:
        print_e("There are issues with the following databases:")
        for i in c:
            print_e(i)
        print_e("Please check your `--database_path` or other supplied database paths.")
        print_e("Try running `install.py` again.")
        sys.exit(1)

# Some functions that print pretty terminal output using ANSI codes
GREEN = "\033[92m"
CYAN = "\033[96m"
WHITE = "\033[97m"
RED = "\033[91m"
BOLD = "\033[1m"
UNDERLINE = "\033[4m"
END = "\033[00m"

def print_tu(text):
    """
    Title,  underlined in green
    """
    print(f"{GREEN} {BOLD} {UNDERLINE} {text} {END}")
def print_t(text):
    """
    Title in green
    """
    print(f"{GREEN} {BOLD} {text} {END}")
def print_h(text):
    """
    Header `[Date Time] text` in cyan
    """
    print(f"{CYAN}{BOLD}[{datetime.datetime.now().strftime('%b %d %I:%M %p')}]{END}",  f"{CYAN} {text} {END}")
def print_n(text):
    """
    Normal `[Date Time] text` in white
    """
    print(f"{WHITE}[{datetime.datetime.now().strftime('%b %d %I:%M %p')}]{END}",  f"{WHITE} {text} {END}")    
def print_e(text):
    """
    Error `[Date Time] text` in red
    """
    print(f"{RED}{BOLD}[{datetime.datetime.now().strftime('%b %d %I:%M %p')}]{END}",  f"{RED} {text} {END}")