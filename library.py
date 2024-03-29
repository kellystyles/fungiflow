import os
import sys
import subprocess
import datetime
import requests
import tarfile
from tqdm import tqdm

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
        if input_args.illumina_f is not None and input_args.illumina_r is not None:
            self.shortf = os.path.abspath(os.path.join(input_args.directory, input_args.illumina_f))
            self.shortr = os.path.abspath(os.path.join(input_args.directory, input_args.illumina_r))
        else:
            self.shortf = None
            self.shortr = None                 
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
                pass
                
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
    print(d)
    if False in d.values():
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
        with open(stdout,  "wt", encoding="utf-8") as out,  open(stderr,  "wt", encoding="utf-8") as err:
            subprocess.run(command, stdout=out, stderr=err, check=False)
    except subprocess.CalledProcessError as error:
        print("An error occurred while running the command:")
        print(f"Command: {error.cmd}")
        print(f"Exit code: {error.returncode}")
        print(f"Output: {error.output}")
        print(f"Error: {error.stderr}")
        print(sys.exc_info()[1])

def execute_shell(command_list, stdout, stderr):
    """"
    Executes a list of bash commands, piping each command in the list after the first etc.
    i.e., command[1] will take the stdout of command[0] as stdin
    
    Commands must be in string format, with spaces included
    """

    # opens the stdout and stderr files for writing
    try:
        with open(stdout,  "wt", encoding="utf-8") as out,  open(stderr,  "wt", encoding="utf-8") as err:
            n = 2
            # If just a single command in the list will execute without PIPE
            if len(command_list) == 1:
                print(command_list[0])
                subprocess.run(command_list[0], stdout=out, stderr=err, shell=True, check=False)
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
    except subprocess.CalledProcessError as error:
        print("An error occurred while running the command:")
        print(f"Command: {error.cmd}")
        print(f"Exit code: {error.returncode}")
        print(f"Output: {error.output}")
        print(f"Error: {error.stderr}")
        print(sys.exc_info()[1])

def download_file(url, save_path):
    """
    Downloads file and displays progress bar. 
    Will check if file is partially downloaded and will resume from last chunk.
    Ensure 'url' and 'save_path' are housed by commas.
    """

    response = None
    total_size = 0
    block_size = 8192  # Default chunk size
    downloaded_size = 0

    if os.path.exists(save_path):
        # Get the size of the existing file
        downloaded_size = os.path.getsize(save_path)
        headers = {"Range": f"bytes={downloaded_size}-"}
        response = requests.get(url, headers=headers, stream=True, timeout=10000)
        response.raise_for_status()
        total_size = int(response.headers.get("content-length")) + downloaded_size
        mode = "ab"  # Append to existing file
    else:
        response = requests.get(url, stream=True,  timeout=10000)
        response.raise_for_status()
        total_size = int(response.headers.get("content-length"))
        mode = "wb"  # Create new file

    with open(save_path, mode) as file, \
        tqdm.tqdm(total=total_size, initial=downloaded_size, unit="B", unit_scale=True, desc=save_path) as progress_bar:
        for chunk in response.iter_content(chunk_size=block_size):
            file.write(chunk)
            progress_bar.update(len(chunk))

    # Check if the downloaded file matches the expected size
    downloaded_size = os.path.getsize(save_path)
    if downloaded_size != total_size:
        raise ValueError("Downloaded file size doesn't match the expected size.")

def extract_tar_gz(file_path, extract_path):
    """
    Extracts a tar.gz file to a given path.
    """
    
    with tarfile.open(file_path, 'r:gz') as tar:
        tar.extractall(extract_path)

def check_databases(input_args, filenames_class_obj):
    """
    Checks if the necessary databases are present based on input arguments.
    Then checks if the database exists and will assign the database to input_args.
    """
    # TRY statements to check what DB paths are supplied
    # In the case where the `--database_path` argument is not supploed, 
    # the supplied individual database paths will be used
    try:
        if input_args.database_path is not None:
            if input_args.its is True:
                its_db = os.path.abspath(os.path.join(input_args.database_path, "ncbi-its"))
            if input_args.kraken2 is True:
                kraken2_db = os.path.abspath(os.path.join(input_args.database_path, "kraken2"))
            if input_args.eggnog is True:
                eggnog_db = os.path.abspath(os.path.join(input_args.database_path, "eggnog"))
        else:
            try:
                its_db = os.path.abspath(input_args.its_db)
            except AttributeError:
                pass
            try: 
                kraken2_db = os.path.abspath(input_args.kraken2_db)
            except AttributeError:
                pass
            try: 
                eggnog_db = os.path.abspath(input_args.eggnog_db)
            except AttributeError:
                pass
    except IndexError:
        print("Please supply a value for `--database_path`")
    except TypeError:
        print("Please supply a value for `--database_path`")
    # Checks if key files are present in each database and assigns the main
    # `input_args` Class the individual databases.
    # Should the database files not exist, as determined by `file_exists_list`,
    # will add this database path to a list.
    c = []
    try:
        ncbi_its = os.path.join(its_db, "ITS_RefSeq_Fungi.nsq")
        taxdb = os.path.join(its_db, "taxdb.bti")
        if file_exists_list([ncbi_its, taxdb], "ITS DB is present", "ITS DB is not present. \
            Please supply a path via `--its_db` or run install.py again.") is False:
            c.append(its_db)
        else:
            filenames_class_obj.its_db = its_db
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
            filenames_class_obj.kraken2_db = kraken2_db    
        print(kraken2_db)
    except UnboundLocalError:
        pass
    try:
        dmnd = os.path.join(eggnog_db, "eggnog_proteins.dmnd")
        taxa = os.path.join(eggnog_db, "eggnog.taxa.db")
        main = os.path.join(eggnog_db, "eggnog.db")
        prots = os.path.join(eggnog_db, "e5.proteomes.faa")
        fungi_dmnd = os.path.join(eggnog_db, "eggnog_prots.fungi.dmnd")
        if file_exists_list([dmnd, taxa, main, prots, fungi_dmnd], "eggnog DB is present", "eggnog DB is not present. \
            Please supply a path via `--eggnog_db` or run install.py again.") is False:
            c.append(eggnog_db)
        else:
            filenames_class_obj.eggnog_db = eggnog_db    
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

# Some functions to manage conda environments for Funannotate
def activate_conda_env(env_name):
    command = f"conda activate {env_name}"
    subprocess.run(command, shell=True, executable="/bin/bash", check=True)

# Function to deactivate the current conda environment
def deactivate_conda_env():
    command = "conda deactivate"
    subprocess.run(command, shell=True, executable="/bin/bash", check=True)

def is_conda_environment_active():
    return "CONDA_PREFIX" in os.environ

def is_conda_environment_available(env_name):
    try:
        # Run 'conda env list' and capture the output
        result = subprocess.run(['conda', 'env', 'list'], capture_output=True, text=True, check=True)
        output_lines = result.stdout.splitlines()

        # Check if the environment name is present in the output
        for line in output_lines:
            if env_name in line:
                return True

        return False
    except subprocess.CalledProcessError:
        # 'conda env list' command returned non-zero exit status
        return False

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