import subprocess
import os
import shlex

class SafeSubprocess:
    @staticmethod
    def run(command, **kwargs):
        # Safely run a command using subprocess
        # Ensure the command is a list and sanitize any shell input
        if isinstance(command, str):
            command = shlex.split(command)
        try:
            result = subprocess.run(command, check=True, capture_output=True, text=True, **kwargs)
            return result.stdout
        except subprocess.CalledProcessError as e:
            print(f"Error occurred: {e}")
            print(f"Command output: {e.output}")
            return None

# Usage:
# output = SafeSubprocess.run('ls -l')
# print(output)