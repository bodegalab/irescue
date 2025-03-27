#!/usr/bin/env python

import argparse
import os
import sys
from multiprocessing import Pool
from functools import partial
from shutil import rmtree
from irescue._version import __version__
from irescue._genomes import __genomes__
from irescue.misc import writerr, versiontuple, run_shell_cmd
from irescue.misc import check_requirement, check_tags
from irescue.map import makeRmsk, getRefs, prepare_whitelist, isec, chrcat
from irescue.map import checkIndex
from irescue.count import split_barcodes, index_features, run_count, formatMM, writeEC

def parseArguments():
    parser = argparse.ArgumentParser(
        prog='IRescue',
        usage="irescue -b <file.bam>"
            " [-g GENOME_ASSEMBLY | -r BED_FILE] [OPTIONS]",
        description="IRescue (Interspersed Repeats single-cell quantifier):"
            " a tool for quantifying transposable elements expression"
            " in scRNA-seq.",
        epilog="Home page: https://github.com/bodegalab/irescue"
    )
    parser.add_argument('-b', '--bam', required=True, metavar='FILE',
                        help="scRNA-seq reads aligned to a reference genome (required).")
    parser.add_argument('-r', '--regions', metavar='FILE',
                        help="Genomic TE coordinates in bed format (at least 4 columns with TE feature name "
                        "(e.g. subfamily) as the 4th column). Takes priority over --genome (default: %(default)s).")
    parser.add_argument('-g', '--genome', metavar='STR', choices=__genomes__.keys(),
                        help="Genome assembly symbol. One of: {} (default: "
                        "%(default)s).".format(', '.join(__genomes__)))
    parser.add_argument('-w', '--whitelist', metavar='FILE',
                        help="Text file of filtered cell barcodes by e.g. Cell Ranger, STARSolo "
                        "or your gene expression quantifier of choice (Recommended. default: %(default)s).")
    parser.add_argument('-c', '--cb-tag', default='CB', metavar='STR',
                        help="BAM tag containing the cell barcode sequence (default: %(default)s).")
    parser.add_argument('-u', '--umi-tag', default='UR', metavar='STR',
                        help="BAM tag containing the UMI sequence (default: %(default)s).")
    parser.add_argument('--no-umi', action='store_true',
                        help="Ignore UMI sequence (for UMI-less technologies, such as SMART-seq).")
    parser.add_argument('-p', '--threads', type=int, default=1, metavar='CPUS',
                        help="Number of cpus to use (default: %(default)s).")
    parser.add_argument('-o', '--outdir', default='irescue_out', metavar='DIR',
                        help="Output directory name (default: %(default)s).")
    parser.add_argument('--min-bp-overlap', type=int, metavar='INT',
                        help="Minimum overlap between read and TE as number of nucleotides (Default: disabled).")
    parser.add_argument('--min-fraction-overlap', type=float, metavar='FLOAT', choices=[x/100 for x in range(101)],
                        help="Minimum overlap between read and TE as a fraction of read's alignment"
                        " (i.e. 0.00 <= NUM <= 1.00) (Default: disabled).")
    parser.add_argument('--max-iters', type=int, metavar='INT', default=100,
                        help="Maximum number of EM iterations (Default: %(default)s).")
    parser.add_argument('--tolerance', type=float, metavar='FLOAT', default=1e-4,
                        help="Log-likelihood change below which convergence is assumed (Default: %(default)s).")
    parser.add_argument('--dump-ec', action='store_true',
                        help="Write a description log file of Equivalence Classes.")
    parser.add_argument('--integers', action='store_true',
                        help="Use if integers count are needed for downstream analysis.")
    parser.add_argument('--samtools', default='samtools', metavar='PATH',
                        help="Path to samtools binary, in case it's not in PATH (Default: %(default)s).")
    parser.add_argument('--bedtools', default='bedtools', metavar='PATH',
                        help="Path to bedtools binary, in case it's not in PATH (Default: %(default)s).")
    parser.add_argument('--no-tags-check', action='store_true',
                        help="Suppress checking for CBtag and UMItag presence in BAM file.")
    parser.add_argument('--keeptmp', action='store_true',
                        help="Keep temporary files under <output_dir>/tmp.")
    parser.add_argument('-v', '--verbose', action='count', default=0,
                        help="Writes additional logging to stderr. Use once for normal verbosity (-v), "
                        "twice for debugging (-vv).")
    parser.add_argument('-V', '--version', action='version', version='%(prog)s {}'.format(__version__),
                        help="Print software's version and exit.")
    return parser


def main():

    # Parse and print arguments
    parser = parseArguments()
    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
    
    if args.no_umi:
        args.umi_tag = ""

    argstr = '\n'.join(f'    {k}: {v}' for k, v in args.__dict__.items())
    sys.stderr.write(f"    IRescue version {__version__}\n{argstr}\n")

    #__tmpdir__ = os.path.join(args.outdir, 'tmp')
    dirs = {
        'out': args.outdir,
        'tmp': os.path.join(args.outdir, 'tmp'),
        'mex': os.path.join(args.outdir, 'counts')
    }


    ####################
    # Preliminar steps #
    ####################

    writerr("Running preliminary checks.")

    # Check requirements
    check_requirement(
        args.bedtools, '2.30.0',
        lambda: versiontuple(
            run_shell_cmd('bedtools --version').split()[1][1:]
        ),
        args.verbose
    )
    check_requirement(
        args.samtools, '1.11',
        lambda: versiontuple(run_shell_cmd('samtools --version').split()[1]),
        args.verbose
    )

    # Check if the selected cell barcode and UMI tags are present in bam file.
    if not args.no_tags_check:
        check_tags(bamFile=args.bam, CBtag=args.cb_tag, UMItag=args.umi_tag,
                   nLines=999999, exit_with_error=True, verbose=args.verbose)

    # Check for bam index file. If not present, will build an index.
    checkIndex(args.bam, verbose=args.verbose)

    # create directories
    for v in dirs.values():
        os.makedirs(v, exist_ok=True)


    ###########
    # Mapping #
    ###########

    writerr("Running mapping step.")

    # set regions object (provided or downloaded bed file)
    regions = makeRmsk(regions=args.regions, genome=args.genome,
                       genomes=__genomes__, tmpdir=dirs['tmp'],
                       outname='rmsk.bed')

    # get list of reference names from bam
    chrNames = getRefs(args.bam, regions)

    # decompress whitelist if compressed
    whitelist = prepare_whitelist(args.whitelist, dirs['tmp'])

    # Allocate threads
    if args.threads > 1:
        pool = Pool(args.threads)

    # Execute intersection between reads and TE coordinates
    writerr(
        "Computing overlap between reads and TEs coordinates in the "
        "following references: {}".format(', '.join(chrNames)),
        level=1, send=args.verbose
    )
    isecFun = partial(
        isec, args.bam, regions, whitelist, args.cb_tag, args.umi_tag,
        args.min_bp_overlap, args.min_fraction_overlap, dirs['tmp'],
        args.samtools, args.bedtools, args.verbose
    )
    if args.threads > 1:
        isecFiles = pool.map(isecFun, chrNames)
    else:
        isecFiles = list(map(isecFun, chrNames))

    # concatenate intersection results
    mappings_file, barcodes_file, features_file = chrcat(
        isecFiles, threads=args.threads, outdir=dirs['mex'],
        tmpdir=dirs['tmp'], bedtools=args.bedtools, verbose=args.verbose
    )


    #########
    # Count #
    #########

    writerr("Running count step.")
    # calculate number of mappings per process
    bc_per_thread = list(split_barcodes(barcodes_file, args.threads))

    # parse features
    feature_index = index_features(features_file)

    # calculate TE counts
    countFun = partial(
        run_count, mappings_file, feature_index, dirs['tmp'], args.no_umi,
        args.dump_ec, args.max_iters, args.tolerance, args.verbose
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
    matrix_files = [ i for i, j in mtxFiles]
    ecdump_files = [ j for i, j in mtxFiles]
    matrix_file = formatMM(
        matrix_files, feature_index, bc_per_thread, dirs['mex']
    )
    writerr(f'Writing sparse matrix to {matrix_file}')
    if args.dump_ec:
        ecdump_file = writeEC(ecdump_files, args.no_umi, outdir=dirs['out'])
        writerr(f'Writing Equivalence Classes to {ecdump_file}')

    if not args.keeptmp:
        writerr('Cleaning up temporary files.', level=1, send=args.verbose)
        rmtree(dirs['tmp'])

    writerr('Done.')
