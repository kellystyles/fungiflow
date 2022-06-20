import os
import subprocess
import datetime

"""
Helper functions for the Fungiflow pipeline.
"""

# Creating the Class for the workflow filepaths
class Files:
    def __init__(self, array):
        self.short_raw_f = os.path.join("raw",f"{array}_R1.fq.gz")
        self.short_raw_r = os.path.join("raw",f"{array}_R2.fq.gz")
    def __str__(self):
        return  str(self.__class__) + '\n' + '\n'.join((str(item) + ' = ' + str(self.__dict__[item]) for item in sorted(self.__dict__)))

def make_path(path):
    """
    Makes a path, regardless of if the parent dir exists or not.
    Functionally equivalent to `mkdir -p` in BASH
    """
    try:
        os.makedirs(path)
    except FileExistsError:
        pass

def file_exists(file,size,success_message,error_message):
    """
    Check if file exists and is not empty.
    Returns error message if the above is not True.
    """
    try:
        if os.stat(file).st_size > int(size):
            print(success_message)
        
    except FileNotFoundError:
        exit(error_message)

def file_exists_bool(file,size,success_message,error_message):
    """
    Check if file exists and is not empty.
    Returns error message if the above is not True.
    """
    try:
        if os.stat(file).st_size > int(size):
            print(success_message)
            return True
    except FileNotFoundError:
        print(error_message)
        return False

def execute(command,stdout,stderr):
    """
    Executes a Bash command and directs the stdout and stderr to separate files, 
    respectively.

    Parameters
    ----------
    command : comma separated list of strings (e.g., ["blastn","-query",infile,"-out",outfile]).
    stdout : file to record standard out.
    stderr : file to record standard error.
    """
    with open(stdout, "wt") as out, open(stderr, "wt") as err:
        subprocess.run(command,stdout=out,stderr=err)

def execute_shell(command,stdout,stderr):
    """
    Executes a Bash command and directs the stdout and stderr to separate files, 
    respectively.

    Parameters
    ----------
    command : comma separated list of strings (e.g., ["blastn","-query",infile,"-out",outfile]).
    stdout : file to record standard out.
    stderr : file to record standard error.
    """
    with open(stdout, "wt") as out, open(stderr, "wt") as err:
        subprocess.run(command,shell=True,stdout=out,stderr=err)
        
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
    Title, underlined in green
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
    print(f"{CYAN}{BOLD}[{datetime.datetime.now().strftime('%b %d %I:%M %p')}]{END}", f"{CYAN} {text} {END}")
def print_n(text):
    """
    Normal `[Date Time] text` in white
    """
    print(f"{WHITE}[{datetime.datetime.now().strftime('%b %d %I:%M %p')}]{END}", f"{WHITE} {text} {END}")    
def print_e(text):
    """
    Error `[Date Time] text` in red
    """
    print(f"{RED}{BOLD}[{datetime.datetime.now().strftime('%b %d %I:%M %p')}]{END}", f"{RED} {text} {END}")