#!/usr/bin/env python

from irescue._version import __version__
from irescue._genomes import __genomes__
from irescue.misc import writerr, check_requirement, versiontuple, run_shell_cmd, check_tags
from irescue.map import makeRmsk, getRefs, prepare_whitelist, isec, chrcat, checkIndex
from irescue.count import split_bc, parse_features, count, formatMM
import argparse, os, sys
from multiprocessing import Pool
from functools import partial
from shutil import rmtree

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
    parser.add_argument('-g', '--genome', default=False, help='Genome assembly symbol. One of: {} (default: False)'.format(', '.join(__genomes__.keys())))
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
    parser.add_argument('--no-tags-check', default=False, action='store_true', help='Suppress checking for CBtag and UMItag presence in bam file (default: False)')
    parser.add_argument('-v', '--verbose', default=False, action='store_true', help='Writes a lot of stuff to stderr, such as chromosomes as they are mapped and cell barcodes as they are processed.')
    parser.add_argument('--version', action='version', version='%(prog)s {}'.format(__version__), help='Print software\'s version and exit')
    return parser

def main():
    
    parser = parseArguments()
    args = parser.parse_args()

    writerr('IRescue job starts')

    # Check requirements
    check_requirement('bedtools', '2.30.0', lambda: versiontuple(run_shell_cmd('bedtools --version').split()[1][1:]), args.verbose)
    check_requirement('samtools', '1.11', lambda: versiontuple(run_shell_cmd('samtools --version').split()[1]), args.verbose)

    # Check if the selected cell barcode and UMI tags are present in bam file.
    if not args.no_tags_check:
        check_tags(bamFile=args.bam, CBtag=args.CBtag, UMItag=args.UMItag,
                   nLines=999999, exit_with_error=True, verbose=args.verbose)

    # Check for bam index file. If not present, will build an index.
    checkIndex(args.bam, verbose=args.verbose)
    
    # create directories
    os.makedirs(args.tmpdir, exist_ok=True)
    os.makedirs(args.outdir, exist_ok=True)

    # set regions object (provided or downloaded bed file)
    regions = makeRmsk(regions=args.regions, genome=args.genome,
                       genomes=__genomes__, tmpdir=args.tmpdir,
                       outname='rmsk.bed')

    # get list of reference names from bam
    chrNames = getRefs(args.bam, regions)

    # decompress whitelist if compressed
    whitelist = prepare_whitelist(args.whitelist, args.tmpdir)

    # Allocate threads
    if args.threads > 1:
        pool = Pool(args.threads)

    # Execute intersection between reads and TE coordinates
    writerr(f"Computing overlap between reads and TEs coordinates in the following references: {', '.join(chrNames)}", send=args.verbose)
    isecFun = partial(
        isec, args.bam, regions, whitelist, args.CBtag, args.UMItag,
        args.tmpdir, args.samtools, args.bedtools, args.verbose
    )
    if args.threads > 1:
        isecFiles = pool.map(isecFun, chrNames)
    else:
        isecFiles = list(map(isecFun, chrNames))

    # concatenate intersection results
    mappings_file, barcodes_file, features_file = chrcat(
        isecFiles, threads=args.threads, outdir=args.outdir,
        tmpdir=args.tmpdir, verbose=args.verbose
    )

    # calculate number of mappings per process
    bc_per_thread = list(split_bc(barcodes_file, args.threads))

    # parse features
    ftlist = dict(parse_features(features_file))

    # calculate TE counts
    countFun = partial(
        count, mappings_file, args.outdir, args.tmpdir, ftlist, args.integers,
        args.verbose
    )
    if args.threads > 1:
        mtxFiles = pool.map(countFun, bc_per_thread)
    else:
        mtxFiles = list(map(countFun, bc_per_thread))
    
    # close processes pool
    if args.threads > 1:
        pool.close()
        pool.join()

    # concatenate matrix files chunks
    matrix_file = formatMM(
        mtxFiles, outdir=args.outdir, features=ftlist, barcodes=bc_per_thread
    )
    writerr(f'Writing sparse matrix to {matrix_file}')

    if not args.keeptmp:
        writerr(f'Cleaning up temporary files.', send=args.verbose)
        rmtree(args.tmpdir)

    writerr('Done.')
