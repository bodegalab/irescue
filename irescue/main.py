#!/usr/bin/env python

from irescue._version import __version__
from irescue._genomes import __genomes__
from irescue.misc import writerr, versiontuple, run_shell_cmd
from irescue.misc import check_requirement, check_arguments, check_tags
from irescue.map import makeRmsk, getRefs, prepare_whitelist, isec, chrcat
from irescue.map import checkIndex
from irescue.count import split_bc, parse_features, count, formatMM, writeEC
import argparse, os
from multiprocessing import Pool
from functools import partial
from shutil import rmtree

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
    parser.add_argument('-b', '--bam',
                        required=True,
                        metavar='FILE',
                        help="scRNA-seq reads aligned to a reference genome "
                        "(required).")
    parser.add_argument('-r', '--regions',
                        metavar='FILE',
                        help="Genomic TE coordinates in bed format. "
                        "Takes priority over --genome (default: %(default)s).")
    parser.add_argument('-g', '--genome',
                        metavar='STR',
                        help="Genome assembly symbol. One of: {} (default: "
                        "%(default)s).".format(', '.join(__genomes__)))
    parser.add_argument('-w', '--whitelist',
                        metavar='FILE',
                        help="Text file of filtered cell barcodes by e.g. "
                        "Cell Ranger, STARSolo or your gene expression "
                        "quantifier of choice (Recommended. "
                        "default: %(default)s).")
    parser.add_argument('-cb', '--CBtag',
                        default='CB',
                        metavar='STR',
                        help="BAM tag containing the cell barcode sequence "
                        "(default: %(default)s).")
    parser.add_argument('-umi', '--UMItag',
                        default='UR',
                        metavar='STR',
                        help="BAM tag containing the UMI sequence "
                        "(default: %(default)s).")
    parser.add_argument('-p', '--threads',
                        type=int,
                        default=1,
                        metavar='CPUS <int>',
                        help="Number of cpus to use (default: %(default)s).")
    parser.add_argument('-o', '--outdir',
                        default='IRescue_out',
                        metavar='DIR',
                        help="Output directory name (default: %(default)s).")
    parser.add_argument('--min-bp-overlap',
                        type=int,
                        metavar='INT',
                        help="Minimum overlap between read and TE as number "
                        "of nucleotides (Default: disabled).")
    parser.add_argument('--min-fraction-overlap',
                        type=float,
                        metavar='FLOAT',
                        help="Minimum overlap between read and TE"
                        " as a fraction of read's alignment"
                        " (i.e. 0.00 <= NUM <= 1.00) (Default: disabled).")
    parser.add_argument('--dumpEC',
                        action='store_true',
                        help="Write a description log file of Equivalence "
                        "Classes.")
    parser.add_argument('--integers',
                        action='store_true',
                        help="Use if integers count are needed for "
                        "downstream analysis.")
    parser.add_argument('--samtools',
                        default='samtools',
                        metavar='PATH',
                        help="Path to samtools binary, in case it's not in "
                        "PATH (Default: %(default)s).")
    parser.add_argument('--bedtools',
                        default='bedtools',
                        metavar='PATH',
                        help="Path to bedtools binary, in case it's not in "
                        "PATH (Default: %(default)s).")
    parser.add_argument('--no-tags-check',
                        action='store_true',
                        help="Suppress checking for CBtag and UMItag "
                        "presence in bam file.")
    parser.add_argument('--keeptmp',
                        action='store_true',
                        help="Keep temporary files.")
    parser.add_argument('--tmpdir',
                        default='IRescue_tmp',
                        metavar='DIR',
                        help="Directory to store temporary files "
                        "(default: %(default)s).")
    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        help="Writes a lot of stuff to stderr, such as "
                        "chromosomes as they are mapped and cell barcodes "
                        "as they are processed.")
    parser.add_argument('-V', '--version',
                        action='version',
                        version='%(prog)s {}'.format(__version__),
                        help="Print software's version and exit.")
    return parser


def main():
    parser = parseArguments()
    args = parser.parse_args()
    args = check_arguments(args)

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
        check_tags(bamFile=args.bam, CBtag=args.CBtag, UMItag=args.UMItag,
                   nLines=999999, exit_with_error=True, verbose=args.verbose)

    # Check for bam index file. If not present, will build an index.
    checkIndex(args.bam, verbose=args.verbose)
    
    writerr('IRescue job starts')
    
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
    writerr(
        "Computing overlap between reads and TEs coordinates in the "
        "following references: {}".format(', '.join(chrNames)),
        send=args.verbose
    )
    isecFun = partial(
        isec, args.bam, regions, whitelist, args.CBtag, args.UMItag,
        args.min_bp_overlap, args.min_fraction_overlap, args.tmpdir,
        args.samtools, args.bedtools, args.verbose
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
        args.dumpEC, args.verbose
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
        matrix_files, outdir=args.outdir, features=ftlist,
        barcodes=bc_per_thread
    )
    writerr(f'Writing sparse matrix to {matrix_file}')
    if args.dumpEC:
        ecdump_file = writeEC(ecdump_files, outdir=args.outdir)
        writerr(f'Writing Equivalence Classes to {ecdump_file}')

    if not args.keeptmp:
        writerr(f'Cleaning up temporary files.', send=args.verbose)
        rmtree(args.tmpdir)

    writerr('Done.')
