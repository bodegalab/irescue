#!/usr/bin/env python

import subprocess
import os
import sys
import gzip
from datetime import datetime
from shutil import which

# Execute a command with subprocess
def run_shell_cmd(cmd):
    p = subprocess.Popen(
        ['/bin/bash', '-o', 'pipefail'],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        preexec_fn=os.setsid,
        text=True
    )
    pid = p.pid
    pgid = os.getpgid(pid)
    stdout, stdin = p.communicate(cmd)
    return stdout.strip('\n')

def check_path(cmdname):
    return which(cmdname) is not None

def versiontuple(v):
    return tuple(map(int, v.split('.')))

# Small function to write a message to stderr with timestamp
def writerr(msg, send=True):
    if send:
        timelog = datetime.now().strftime("%m/%d/%Y - %H:%M:%S")
        message = f'[{timelog}] '
        if not msg[-1]=='\n':
            msg += '\n'
        message += msg
        sys.stderr.write(message)

# Test if file is gzip
def testGz(input_file):
    with gzip.open(input_file, 'rb') as f:
        try:
            f.read(1)
            return True
        except gzip.BadGzipFile:
            return False

# Uncompress gzipped file
def unGzip(input_file, output_file):
    with gzip.open(input_file, 'rb') as fin,\
    open(output_file, 'w') as fout:
        for line in fin:
            fout.write(line.decode('utf-8'))
    return output_file

# Count number of lines in file (also compressed)
def getlen(file):
    if testGz(file):
        f = gzip.open(file, 'rb')
    else:
        f = open(file, 'r')
    out = sum(1 for line in f)
    f.close()
    return out

# Flatten a list of sublists
def flatten(x):
    return [item for sublist in x for item in sublist]
