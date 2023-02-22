#!/usr/bin/env python

import subprocess
import os
import sys
import gzip
from datetime import datetime
from shutil import which
import pysam

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

def check_requirement(cmd, required_version, parser, verbose):
    if not check_path(cmd):
        writerr(f"ERROR: Couldn't find {cmd} in PATH. Please install {cmd} >={required_version} and try again.", error=True)
    else:
        try:
            version = parser()
            if version < versiontuple(required_version):
                writerr(f"WARNING: Found {cmd} version {version}. Versions prior {required_version} are not supported.")
            else:
                writerr(f"Found {cmd} version {version}. Proceeding.", send=verbose)
        except:
            writerr(f"WARNING: Found {cmd} but couldn't parse its version. NB: {cmd} versions prior {required_version} are not supported.")

# Small function to write a message to stderr with timestamp
def writerr(msg, error=False, send=True):
    if send:
        timelog = datetime.now().strftime("%m/%d/%Y - %H:%M:%S")
        message = f'[{timelog}] '
        if not msg[-1]=='\n':
            msg += '\n'
        message += msg
        if error:
            sys.exit(message)
        else:
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

# Check if bam contains barcode and umi tags
def check_tags(bamFile, CBtag, UMItag, nLines=False, exit_with_error=True, verbose=False):
    with pysam.AlignmentFile(bamFile, 'rb') as f:
        c = 1
        writerr(f"Testing bam file for {CBtag} and {UMItag} tags presence. Will stop at the first occurrence.", send=verbose)
        for read in f:
            if nLines and c >= nLines:
                break
            elif c % 1000000 == 0:
                writerr(f"WARNING: Couldn't find {CBtag} and {UMItag} tags in the first {c} records. Did you select the right tags? Continuing parsing bam...")
            try:
                read.get_tag(CBtag) and read.get_tag(UMItag)
                writerr(f"Found {CBtag} and {UMItag} tags occurrence in bam's line {c}.")
                return(True)
            except:
                c += 1
                pass
    if exit_with_error:
        writerr(
            """
            ERROR: Couldn't find {} and {} tags in {}the bam file.
            Check you bam file for the presence of tags for cell barcode
            and UMI sequences, then provide them to IRescue through
            --CBtag and --UMItag flags.
            If you expect few alignments to contain the tags, you can
            suppress this check with --no-tags-check
            """.format(
                CBtag,
                UMItag,
                f'the first {nLines} lines of ' if nLines else ''
            ),
            error=True
        )
    else:
        return(False)
