import os
import sys
import subprocess
import datetime

"""
Helper functions for the Fungiflow pipeline.
"""

# Creating the Class for the workflow filepaths
class Files:
    def __init__(self, input_args):
        # defines paths for short reads
        self.shortf = os.path.abspath(os.path.join(input_args.directory, input_args.illumina_f))
        self.shortr = os.path.abspath(os.path.join(input_args.directory, input_args.illumina_r))     
        # defines whether singularity commands are needed to execute the pipeline
        if input_args.singularity_image is not None: 
            self.singularity = ["singularity", "exec", input_args.singularity_image]
        else:
            self.singularity = ""
    def __str__(self):
        return  str(self.__class__) + '\n' + '\n'.join((str(item) + ' = ' + str(self.__dict__[item]) for item in sorted(self.__dict__)))
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
        print_n(success_message)
        return True
    else:
        if len(error_message) > 1:
            print_e(error_message)
        return False

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
    print(" ".join(command))
    with open(stdout,  "wt") as out,  open(stderr,  "wt") as err:
        subprocess.run(command, stdout=out, stderr=err)

def execute_shell(command_list, stdout, stderr):
    """"
    Executes a list of bash commands, piping each command in the list after the first etc.
    i.e., command[1] will take the stdout of command[0] as stdin
    
    Commands must be in string format, with spaces included
    """
    
    # opens the stdout and stderr files for writing
    with open(stdout,  "wt") as out,  open(stderr,  "wt") as err:
        n = 2
        # takes first command in command_list as first command, notice no `input=`
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