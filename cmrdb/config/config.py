import json
import os
import sys
import subprocess

"""
Set the default conda environment path.
"""

def run(command, stdin=None):
    '''
    Run a subprocess.check_output() with the given command with
    'bash -c command'
    returning the stdout. If the command fails (i.e. has a non-zero exitstatus),
    raise a ExternCalledProcessError that includes the $stderr as part of
    the error message
    Parameters
    ----------
    command: str
        command to run
    stdin: str or None
        stdin to be provided to the process, to subprocess.communicate.
    Returns
    -------
    Standard output of the run command
    Exceptions
    ----------
    extern.ExternCalledProcessError including stdout and stderr of the run
    command should it return with non-zero exit status.
    '''
    #logging.debug("Running extern cmd: %s" % command)

    using_stdin = stdin is not None
    process = process = subprocess.Popen(
        ["bash",'-o','pipefail',"-c", command],
        stdin= (subprocess.PIPE if using_stdin else None),
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate(stdin)

    if process.returncode != 0:
        raise ExternCalledProcessError(process, command, stdout.decode(), stderr.decode())
    return stdout

def get_conda_path():
    try:
        CONDA_PATH = run("conda config --show envs_dirs | sed -n 2p | sed 's/  - //;s/$/\//'").decode().rstrip()
        return CONDA_PATH
    except KeyError:
        print('\n' + '=' * 80)
        print(' ERROR '.center(80))
        print('_' * 80 + '\n')
        print("The 'CONDA_ENV_PATH' environment variable is not defined.".center(80) + '\n')
        print('Please set this variable to your default server/home directory conda environment path.'.center(80))
        print('Alternatively, use --conda-prefix flag.'.center(80))
        print('=' * 80)
        sys.exit(1)
