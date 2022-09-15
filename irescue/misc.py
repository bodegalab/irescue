#!/usr/bin/env python

import subprocess
import os
from datetime import datetime
import sys
import gzip

def run_shell_cmd(cmd):
    '''Wrapper function to execute subprocesses'''
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

def time():
    '''Small function to write out the current date and time'''
    return datetime.now().strftime("%m/%d/%Y - %H:%M:%S")

# Small function to write a message to stderr with timestamp
def writerr(msg, send=True):
    if send:
        message = f'[{time()}] '
        if not msg[-1]=='\n':
            msg += '\n'
        message += msg
        sys.stderr.write(message)

def testGz(input_file):
    '''Test if file is gzip'''
    with gzip.open(input_file, 'rb') as f:
        try:
            f.read(1)
            return True
        except gzip.BadGzipFile:
            return False

def unGzip(input_file, output_file):
    '''Decompress gzip files'''
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
    '''Flattens a list of sublists'''
    return [item for sublist in x for item in sublist]
