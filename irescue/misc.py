#!/usr/bin/env python

import subprocess
import os
import sys
import gzip
from datetime import datetime
from shutil import which
import pysam

def run_shell_cmd(cmd):
    """
    Execute a command on bash shell with subprocess.
    """
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
    """
    Check if a command name is available in PATH.
    """
    return which(cmdname) is not None

def versiontuple(version):
    """
    Convert a semver string "X.Y.Z" to a tuple of integers (X, Y, Z).
    """
    return tuple(map(int, version.split('.')))

def check_requirement(cmd, required_version, parser, verbose):
    """
    Check if the required version for a software has been installed.
    """
    if not check_path(cmd):
        writerr(
            f"ERROR: Couldn't find {cmd} in PATH. Please install "
            f"{cmd} >={required_version} and try again.",
            error=True
        )
    else:
        try:
            version = parser()
            if version < versiontuple(required_version):
                writerr(
                    f"WARNING: Found {cmd} version {version}. "
                    f"Versions prior {required_version} are not supported."
                )
            else:
                writerr(
                    f"Found {cmd} version {version}. Proceeding.",
                    level=1, send=verbose
                )
        except:
            writerr(
                f"WARNING: Found {cmd} but couldn't parse its version. "
                f"NB: {cmd} versions prior {required_version} are "
                f"not supported."
            )

# Small function to write a message to stderr with timestamp
def writerr(msg, error=False, level=0, send=0):
    """
    Write a message to stderr with timestamp.

    Parameters
    ----------
    msg: string
        Message to write to stderr.
    error: bool
        Set True if the message is an error (write with sys.exit).
    level: int
        Verbosity level of the message.
    send: int
        If >=verbosity, the message will be sent.
    """
    if send>=level or error:
        timelog = datetime.now().strftime("%Y/%m/%d - %H:%M:%S")
        message = f'[{timelog}] '
        if not msg[-1]=='\n':
            msg += '\n'
        message += msg
        if error:
            sys.exit(message)
        else:
            sys.stderr.write(message)

def testGz(input_file):
    """
    Check if a file is gzip compressed.
    """
    with gzip.open(input_file, 'rb') as f:
        try:
            f.read(1)
            return True
        except gzip.BadGzipFile:
            return False

# Uncompress gzipped file
def unGzip(input_file, output_file):
    """
    Uncompress a gzip file to a new file.

    Parameters
    ----------
    input_file: gzip file to decompress.
    output_file: uncompressed file to write.
    """
    with gzip.open(input_file, 'rb') as fin,\
    open(output_file, 'w') as fout:
        for line in fin:
            fout.write(line.decode('utf-8'))
    return output_file

def getlen(file):
    """
    Count the number of lines in a file (plain or gzip-compressed).
    """
    if testGz(file):
        f = gzip.open(file, 'rb')
    else:
        f = open(file, 'r')
    out = sum(1 for line in f)
    f.close()
    return out

def check_tags(
        bamFile, CBtag, UMItag,
        nLines=None, exit_with_error=True, verbose=False
):
    """
    Check if BAM file contains barcode and UMI tags.
    Stops at first co-presence of both tags.

    Parameters
    ----------
    bamFile: BAM file.
    CBtag: string
        Cell Barcode sequence SAM tag
    UMItag: string
        UMI sequence SAM tag
    nLines: int
        Check first N lines.
    exit_with_error: bool
        Return a sys.exit message if tags are not found.
    verbose: bool
        Write progress info to stderr.
    """
    with pysam.AlignmentFile(bamFile, 'rb') as f:
        c = 1
        if not UMItag:
            writerr(
                f"Testing BAM file for {CBtag} tag presence. "
                "Will stop at the first occurrence.",
                level=1, send=verbose
            )
        else:
            writerr(
                f"Testing BAM file for {CBtag} and {UMItag} tags presence. "
                "Will stop at the first occurrence.",
                level=1, send=verbose
            )    
        for read in f:
            if nLines and c >= nLines:
                break
            elif c % 1000000 == 0:
                if not UMItag:     ## could have more sense to have two for loops for the two conditions (so to not check at each iteration but just once)
                    writerr(
                        f"WARNING: Couldn't find {CBtag} tag in "
                        f"the first {c} records. Did you select the right tag? "
                        "Continuing parsing bam until first occurrence..."
                    )
                else:
                    writerr(
                        f"WARNING: Couldn't find {CBtag} and {UMItag} tags in "
                        f"the first {c} records. Did you select the right tags? "
                        "Continuing parsing bam until first occurrence..."
                    )    
            try:
                if not UMItag:
                    read.get_tag(CBtag)
                    writerr(
                        f"Found {CBtag} tag occurrence "
                        f"in BAM's line {c}."
                    )
                else:
                    read.get_tag(CBtag) and read.get_tag(UMItag)
                    writerr(
                        f"Found {CBtag} and {UMItag} tags occurrence "
                        f"in BAM's line {c}."
                    )
                return(True)
            except:
                c += 1
                pass
    if exit_with_error:
        if not UMItag:
            writerr(
                """
                ERROR: Couldn't find {} tag in {}the BAM file.
                Check your BAM file for the presence of tag for cell barcode,
                 then provide it to IRescue through --CBtag flag.
                If you expect few alignments to contain the tag, you can
                suppress this check with --no-tags-check
                """.format(
                    CBtag,
                    f'the first {nLines} lines of ' if nLines else ''
                ),
                error=True
            )
        else:
            writerr(
                """
                ERROR: Couldn't find {} and {} tags in {}the BAM file.
                Check your BAM file for the presence of tags for cell barcode
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

def get_ranges(num, div):
    """
    Splits an integer X into N integers whose sum is equal to X.
    """
    split = int(num/div)
    for i in range(0, num, split):
        j = i + split
        if j > num-split:
            j = num
            yield range(i, j)
            break
        yield range(i, j)
