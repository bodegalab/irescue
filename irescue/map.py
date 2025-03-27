#!/usr/bin/env python

import requests
import io
import os
from pysam import idxstats, AlignmentFile, index
from gzip import open as gzopen
from irescue.misc import testGz
from irescue.misc import writerr
from irescue.misc import unGzip
from irescue.misc import run_shell_cmd
from irescue.misc import getlen

# Check if bam file is indexed
def checkIndex(bamFile, verbose):
    with AlignmentFile(bamFile) as bam:
        if not bam.has_index():
            writerr('BAM index not found. Attempting to index the BAM...')
            try:
                index(bamFile)
            except Exception as e:
                writerr(
                    "ERROR: Couldn't index the BAM file. Is your BAM file "
                    "sorted? If not, please sort it by coordinate "
                    f"with `samtools index {bamFile}`.\n\n{e}",
                    error=True
                )
            else:
                writerr('BAM indexing done.')
        else:
            if verbose:
                writerr(f'Found index for BAM file {bamFile}.',
                        level=1, send=verbose)

# Check repeatmasker regions bed file format. Download if not provided.
# Returns the path of the repeatmasker bed file.
def makeRmsk(regions, genome, genomes, tmpdir, outname):
    # if a repeatmasker bed file is provided, use that
    if regions:
        if testGz(regions):
            f = gzopen(regions, 'rb')
            def rl(x):
                return x.readline().decode()
        else:
            f = open(regions, 'r')
            def rl(x):
                return x.readline()
        # skip header
        line = rl(f)
        while line[0] == '#':
            line = rl(f)
        # check for minimum column number
        if len(line.strip().split('\t')) < 4:
            writerr(
                "Error: please provide a tab-separated BED file with at "
                "least 4 columns and TE feature name (e.g. subfamily) "
                "in 4th column.",
                error=True
            )
        f.close()
        out = regions
    # if no repeatmasker file is provided, and a genome assembly name is
    # provided, download and prepare a rmsk.bed file
    elif genome:
        url, header_lines = genomes[genome]
        writerr(
            "Downloading and parsing RepeatMasker annotation for "
            f"assembly {genome} from {url} ...")
        try:
            response = requests.get(url, stream=True, timeout=60)
        except Exception as e:
            writerr(
                "ERROR: Download of RepeatMasker annotation failed. "
                f"Couldn't connect to host.\n\n{e}",
                error=True
            )
        rmsk = gzopen(io.BytesIO(response.content), 'rb')
        out = os.path.join(tmpdir, outname)
        with open(out, 'w') as f:
            # print header
            h = ['#chr','start','end','name','score','strand']
            h = '\t'.join(h)
            h += '\n'
            f.write(h)
            # skip rmsk header
            for _ in range(header_lines):
                next(rmsk)
            # parse rmsk
            fams_to_skip = [
                'Low_complexity',
                'Simple_repeat',
                'rRNA',
                'scRNA',
                'srpRNA',
                'tRNA'
            ]
            for line in rmsk:
                lst = line.decode('utf-8').strip().split()
                strand, subfamily, famclass = lst[8:11]
                if famclass.split('/')[0] in fams_to_skip:
                    continue
                # concatenate family and class with subfamily
                subfamily += '#' + famclass
                score = lst[0]
                chr, start, end = lst[4:7]
                # make coordinates 0-based
                start = str(int(start)-1)
                if strand != '+':
                    strand = '-'
                outl = '\t'.join([chr, start, end, subfamily, score, strand])
                outl += '\n'
                f.write(outl)
    else:
        writerr(
            "Error: it is mandatory to define either --regions OR "
            "--genome parameter.",
            error=True
        )
    return(out)

# Uncompress the whitelist file if compressed.
# Return the whitelist path, or False if not using a whitelist.
def prepare_whitelist(whitelist, tmpdir):
    if whitelist and testGz(whitelist):
        wlout = os.path.join(tmpdir, 'whitelist.tsv')
        whitelist = unGzip(whitelist, wlout)
    return whitelist

# Get list of reference names from BAM file, skipping those without reads.
def getRefs(bamFile, bedFile):
    chrNames = list()
    for line in idxstats(bamFile).strip().split('\n'):
        fields = line.strip().split('\t')
        if int(fields[2])>0:
            chrNames.append(fields[0])
    bedChrNames = set()
    if testGz(bedFile):
        with gzopen(bedFile, 'rb') as f:
            for line in f:
                bedChrNames.add(line.decode().split('\t')[0])
    else:
        with open(bedFile, 'r') as f:
            for line in f:
                bedChrNames.add(line.split('\t')[0])
    skipChr = [x for x in chrNames if x not in bedChrNames]
    if skipChr:
        writerr(
            "WARNING: The following references contain read alignments but "
            "are not found in the TE annotation and will be skipped: "
            f"{', '.join(skipChr)}")
        chrNames = [x for x in chrNames if x in bedChrNames]
    if chrNames:
        return chrNames
    else:
        writerr(
            """
            ERROR: Reference names not matching between BAM and TE annotation.
            If your BAM follows the ENSEMBL nomenclature (i.e. 1, 2, etc...),
            you can either change it to UCSC (chr1, chr2, etc...), or use a
            custom TE annotation with ENSEMBL chromosome names.
            """,
            error=True
        )

# Intersect reads with repeatmasker regions. Return the intersection file path.
def isec(bamFile, bedFile, whitelist, CBtag, UMItag, bpOverlap, fracOverlap,
         tmpdir, samtools, bedtools, verbose, chrom):
    refdir = os.path.join(tmpdir, 'refs')
    isecdir = os.path.join(tmpdir, 'isec')
    os.makedirs(refdir, exist_ok=True)
    os.makedirs(isecdir, exist_ok=True)

    refFile = os.path.join(refdir, chrom + '.bed.gz')
    isecFile = os.path.join(isecdir, chrom + '.isec.txt.gz')

    # split bed file by chromosome
    sort = 'LC_ALL=C sort -k1,1 -k2,2n --buffer-size=1G'
    if bedFile[-3:] == '.gz':
        cmd0 = f'zcat {bedFile} | gawk \'$1=="{chrom}"\' '
        cmd0 += f' | {sort} | gzip > {refFile}'
    else:
        cmd0 = f'gawk \'$1=="{chrom}"\' {bedFile} | {sort} | gzip > {refFile}'

    # command streaming alignments for intersection
    if whitelist:
        stream = f' <({samtools} view -h {bamFile} -D {CBtag}:{whitelist} {chrom} | '
    else:
        stream = f' <({samtools} view -h {bamFile} {chrom} | '
    stream += ' gawk \'!($1~/^@/) { split("", tags); '
    stream += ' for (i=12;i<=NF;i++) {split($i,tag,":"); tags[tag[1]]=tag[3]}; '
    # Discard records without CB tag, unvalid STARSolo CBs, missing UMI tag,
    # UMIs with Ns and homopolymer UMIs
    if UMItag:
        stream += f' if(tags["{CBtag}"]~/^(|-)$/ || tags["{UMItag}"]~/.*N.*/ || '
        stream += f' tags["{UMItag}"]~/^$|^(A+|G+|T+|C+)$/) {{next}}; '
        # Append CB and UMI to read name
        stream += f' $1=$1"/"tags["{CBtag}"]"/"tags["{UMItag}"]; '
    else:
        stream += f' if(tags["{CBtag}"]~/^(|-)$/) {{next}}; '
        # Append CB to read name, and read name again (replacing UMI)
        stream += f' $1=$1"/"tags["{CBtag}"]"/"$1; '
    stream += ' } '
    stream += ' { OFS="\\t"; print }\' | '
    stream += f' {samtools} view -u - | '
    stream += f' {bedtools} bamtobed -i stdin -bed12 -split -splitD) '

    # filter by minimum overlap between read and feature, if set
    ovfrac = f' -f {fracOverlap} ' if fracOverlap else ''
    ovbp = f' $NF>={bpOverlap} ' if bpOverlap else ''

    # intersection command
    cmd = f'{bedtools} intersect -a {stream} -b {refFile} '
    cmd += f' -split -bed -wo -sorted {ovfrac} | gawk -vOFS="\\t" \'{ovbp} '
    # remove mate information from read name
    cmd += ' { sub(/\\/[12]$/,"",$4); '
    # concatenate CB and UMI with feature name
    cmd += ' n=split($4,qname,/\\//); '
    cmd += ' print qname[n-1]"\\t"qname[n]"\\t"qname[1]"\\t"$16 }\' '
    cmd += f' | gzip > {isecFile}'

    writerr(f'Extracting {chrom} reference', level=2, send=verbose)
    run_shell_cmd(cmd0)
    writerr(f'Mapping alignments to {chrom}', level=1, send=verbose)
    run_shell_cmd(cmd)
    writerr(f'Mapped {chrom}', level=1, send=verbose)

    return isecFile

# Concatenate and sort data obtained from isec()
def chrcat(filesList, threads, outdir, tmpdir, bedtools, verbose):
    os.makedirs(outdir, exist_ok=True)
    mappings_file = os.path.join(tmpdir, 'mappings.tsv.gz')
    barcodes_file = os.path.join(outdir, 'barcodes.tsv.gz')
    features_file = os.path.join(outdir, 'features.tsv.gz')
    bedFiles = ' '.join(filesList)
    sort_threads = int(threads / 2 - 1)
    sort_threads = sort_threads if sort_threads>0 else 1

    # sort and summarize UMI-READ-TE mappings
    # (if --no-umi, UMI is replaced by READ)
    sort_res = f'--parallel {sort_threads} --buffer-size 2G'
    cmd0 = f'zcat {bedFiles}'
        # input: "CB UMI READ FEAT"
    cmd0 += f' | LC_ALL=C sort -u {sort_res}'
    cmd0 += f' | {bedtools} groupby -g 1,2,3 -c 4 -o distinct'
        # result: "CB UMI READ FEATs"
    cmd0 += f' | LC_ALL=C sort -k1,2 -k4,4 {sort_res}'
    cmd0 += f' | {bedtools} groupby -g 1,2,4 -c 3 -o count_distinct'
        # result: "CB UMI FEATs count"
    cmd0 += f' | gzip > {mappings_file}'

    # write barcodes.tsv.gz file
    cmd1 = f'zcat {mappings_file} | cut -f1 | uniq | gzip > {barcodes_file} '

    # write features.tsv.gz file
    cmd2 = f'zcat {mappings_file} '
    cmd2 += ' | cut -f3 | sed \'s/,/\\n/g\' | gawk \'!x[$1]++ { '
    cmd2 += ' print $1"\\t"gensub(/#.+/,"",1,$1)"\\tGene Expression" }\' '
    cmd2 += f' | LC_ALL=C sort -u | gzip > {features_file} '

    writerr('Concatenating mappings', level=1, send=verbose)
    run_shell_cmd(cmd0)
    if getlen(mappings_file) == 0:
        writerr(
            f'No read-TE mappings found in {mappings_file}.'
            ' Check annotation and temporary files to troubleshoot.',
            error=True
        )
    writerr(f'Writing mapped barcodes to {barcodes_file}')
    run_shell_cmd(cmd1)
    if getlen(barcodes_file) == 0:
        writerr(
            f'No features written in {features_file}.'
            ' Check BAM format and reference annotation (e.g. chr names)'
            ' to troubleshoot.',
            error=True
        )
    writerr(f'Writing mapped features to {features_file}')
    run_shell_cmd(cmd2)
    if getlen(features_file) == 0:
        writerr(
            f'No features written in {features_file}.'
            ' Check annotation and temporary files to troubleshoot.',
            error=True
        )

    return mappings_file, barcodes_file, features_file