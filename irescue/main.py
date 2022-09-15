#!/usr/bin/env python

import argparse
from shutil import rmtree
from importlib.metadata import version
from irescue._genomes import __genomes__
from irescue.misc import writerr
from irescue.map import makeRmsk, getRefs, prepare_whitelist, isec, chrcat, checkIndex
from irescue.count import split_bc, parse_features, count, formatMM
import multiprocessing
import os
from functools import partial

def parseArguments():
    parser = argparse.ArgumentParser(
        prog='IRescue',
        usage='''
Automatically download TE annotation:
    irescue -b <file.bam> -g <genome_assembly> [OPTIONS]
Provide a custom TE annotation:
    irescue -b <file.bam> -r <repeatmasker.bed> [OPTIONS]
''',
        description='''IRescue (Interspersed Repeats single-cell quantifier):
a tool for quantifying tansposable elements expression in scRNA-seq.
''',
        epilog='Home page: https://github.com/bodegalab/irescue'
    )
    parser.add_argument('-b','--bam', required=True, help='scRNA-seq reads aligned to a reference genome')
    parser.add_argument('-r', '--regions', default=False, help='Genomic TE coordinates in bed format. Takes priority over --genome paramter (default: False).')
    parser.add_argument('-g', '--genome', default=False, help='Genome assembly symbol. One of: {} (default: False)'.format(','.join(__genomes__.keys())))
    parser.add_argument('-p','--threads', type=int, default=1, help='Number of cpus to use (default: 1)')
    parser.add_argument('-w','--whitelist', default=False, help='Text file of filtered cell barcodes, e.g. by Cell Ranger, STARSolo or your gene expression quantifier of choice (Recommended. Default: False)')
    parser.add_argument('--CBtag', type=str, default='CB', help='BAM tag containing the cell barcode sequence (default: CB)')
    parser.add_argument('--UMItag', type=str, default='UR', help='BAM tag containing the UMI sequence (default: UR)')
    parser.add_argument('--integers', default=False, action='store_true', help='Use if integers count are needed for downstream analysis (default: False)')
    parser.add_argument('--outdir', type=str, default='./IRescue_out/', help='Output directory name (default: IRescue_out)')
    parser.add_argument('--tmpdir', type=str, default='./IRescue_tmp/', help='Directory to store temporary files (default: IRescue_tmp)')
    parser.add_argument('--keeptmp', default=False, action='store_true', help='Keep temporary files (default: False).')
    parser.add_argument('--samtools', type=str, default='samtools', help='Path to samtools binary, in case it\'s not in PATH (Default: samtools)')
    parser.add_argument('--bedtools', type=str, default='bedtools', help='Path to bedtools binary, in case it\'s not in PATH (Default: bedtools)')
    parser.add_argument('-v', '--verbose', default=False, action='store_true', help='Writes a lot of stuff to stderr, such as chromosomes as they are mapped and cell barcodes as they are processed.')
    parser.add_argument('--version', action='version', version='%(prog)s {}'.format(version(parser.prog)), help='Print software\'s version and exit')
    return parser

def main():
    
    parser = parseArguments()
    args = parser.parse_args()

    writerr('IRescue job starts')

    # create directories (TODO: make function to take care of all dirs)
    os.makedirs(args.tmpdir, exist_ok=True)
    os.makedirs(args.outdir, exist_ok=True)

    # check if bam file is indexed. If not, build an index.
    checkIndex(args.bam, verbose=args.verbose)

    # set regions object (provided or downloaded bed file)
    regions = makeRmsk(regions=args.regions,
                       genome=args.genome,
                       genomes=__genomes__,
                       tmpdir=args.tmpdir,
                       outname='rmsk.bed')

    # get list of reference names from bam
    chrNames = getRefs(args.bam)

    # decompress whitelist if compressed
    whitelist = prepare_whitelist(args.whitelist, args.tmpdir)

    # multiprocess execution of isec()
    pool = multiprocessing.Pool(args.threads)
    func = partial(
        isec, args.bam, regions, whitelist, args.CBtag, args.UMItag,
        args.tmpdir, args.samtools, args.bedtools, args.verbose
    )
    isecFiles = pool.map(func, chrNames)
    # concatenate results
    mappings_file, barcodes_file, features_file = chrcat(
        isecFiles, threads=args.threads, outdir=args.outdir,
        tmpdir=args.tmpdir, verbose=args.verbose
    )

    # calculate number of mappings per process
    bc_per_thread = list(split_bc(barcodes_file, args.threads))

    # parse features
    ftlist = dict(parse_features(features_file))

    # calculate TE counts
    #counts = count(mappings_file, outdir=args.outdir, intcount=args.integers)
    countFun = partial(
        count, mappings_file, args.outdir, args.tmpdir, ftlist, args.integers, args.verbose
    )
    mtxFiles = pool.map(countFun, bc_per_thread)

    # concatenate matrix files chunks
    matrix_file = formatMM(
        mtxFiles, outdir=args.outdir, features=ftlist, barcodes=bc_per_thread
    )
    writerr(f'Writing sparse matrix to {matrix_file}')

    if not args.keeptmp:
        writerr(f'Cleaning up temporary files.', args.verbose)
        rmtree(args.tmpdir)

    writerr('Done.')
