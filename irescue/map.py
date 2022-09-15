#!/usr/bin/env python

from irescue.misc import testGz
from irescue.misc import writerr
from irescue.misc import unGzip
from irescue.misc import run_shell_cmd
from pysam import idxstats, AlignmentFile, index
from gzip import open as gzopen
from os import makedirs
import sys, requests, io

# Check if bam file is indexed
def checkIndex(bamFile, verbose):
    with AlignmentFile(bamFile) as bam:
        if not bam.has_index():
            writerr('BAM index not found. Attempting to index the BAM...')
            try:
                index(bamFile)
            except:
                sys.exit(f'Couldn\'t index the BAM file. Please do so manually with `samtools index {bamFile}`.')
            else:
                writerr('BAM indexing done.')
        else:
            if verbose:
                writerr(f'Found index for BAM file {bamFile}', verbose)


# Check repeatmasker regions bed file format. Download if not provided.
# Returns the path of the repeatmasker bed file.
def makeRmsk(regions, genome, genomes, tmpdir, outname):
    # if a repeatmasker bed file is provided, use that
    if regions:
        if testGz(regions):
            f = gzopen(regions, 'rb')
            rl = lambda x: x.readline().decode()
        else:
            f = open(regions, 'r')
            rl = lambda x: x.readline()
        # skip header
        line = rl(f)
        while line[0] == '#':
            line = rl(f)
        # check for minimum column number
        if len(line.strip().split('\t')) < 4:
            sys.exit('Error: please provide a tab-separated BED file with at least 4 columns and TE feature name (e.g. subfamily) in 4th column.')
        f.close()
        out = regions
    # if no repeatmasker file is provided, and a genome assembly name is provided, download and prepare a rmsk.bed file
    elif genome:
        url, header_lines = genomes[genome]
        writerr(f'Downloading and parsing RepeatMasker annotation for assembly {genome} from {url}...')
        response = requests.get(url)
        rmsk = gzopen(io.BytesIO(response.content), 'rb')
        out = tmpdir + '/' + outname
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
            for line in rmsk:
                lst = line.decode('utf-8').strip().split()
                strand, subfamily, famclass = lst[8:11]
                if famclass.split('/')[0] in ['Low_complexity','Simple_repeat','rRNA','scRNA','srpRNA','tRNA']:
                    continue
                # concatenate family and class with subfamily
                subfamily += '~' + famclass
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
        sys.exit('Error: it is mandatory to define either --regions OR --genome paramter.')
    return(out)

# Uncompress the whitelist file if compressed.
# Return the whitelist path, or False if not using a whitelist.
def prepare_whitelist(whitelist, tmpdir):
    if whitelist and testGz(whitelist):
        wlout = tmpdir + '/whitelist.tsv'
        whitelist = unGzip(whitelist, wlout)
    return whitelist

# Get list of reference names from BAM file, skipping those without reads.
def getRefs(bamFile):
    chrNames = list()
    for line in idxstats(bamFile).strip().split('\n'):
        l = line.strip().split('\t')
        if int(l[2])>0:
            chrNames.append(l[0])
    return chrNames

# Intersect reads with repeatmasker regions. Return the intersection file path.
def isec(bamFile, bedFile, whitelist, CBtag, UMItag, tmpdir, samtools, bedtools, verbose, chrom):
    refdir = tmpdir + '/refs/'
    isecdir = tmpdir + '/isec/'
    makedirs(refdir, exist_ok=True)
    makedirs(isecdir, exist_ok=True)

    refFile = refdir + chrom + '.bed.gz'
    isecFile = isecdir + chrom + '.isec.bed.gz'

    # split bed file by chromosome
    sort = 'LC_ALL=C sort -k1,1 -k2,2n --buffer-size=1G'
    if bedFile[-3:] == '.gz':
        cmd0 = f'zcat {bedFile} | awk \'$1=="{chrom}"\' | {sort} | gzip > {refFile}'
    else:
        cmd0 = f'awk \'$1=="{chrom}"\' {bedFile} | {sort} | gzip > {refFile}'

    # command streaming alignments for intersection
    if whitelist:
        stream = f' <({samtools} view -h {bamFile} -D {CBtag}:{whitelist} {chrom} | '
    else:
        stream = f' <({samtools} view -h {bamFile} {chrom} | '
    stream += ' awk \'!($1~/^@/) {'
    stream += ' for (i=12;i<=NF;i++) {split($i,tag,":"); tags[tag[1]]=tag[3]}; '
    # discard unvalid STARSolo CBs and homopolymer UMIs
    stream += f' if(tags["{CBtag}"]=="-" || tags["{UMItag}"]~/^(A+|G+|T+|C+|N+)$/) {{next}}; '
    # append CB and UMI to read name
    stream += f' $1=$1"/"tags["{CBtag}"]"/"tags["{UMItag}"]; '
    stream += ' } '
    stream += ' { OFS="\\t"; print }\' | '
    stream += f' {samtools} view -u - | '
    stream += f' {bedtools} bamtobed -i stdin -bed12 -split -splitD) '

    # intersection command
    cmd = f'{bedtools} intersect -a {stream} -b {refFile} -split -bed -wo -sorted | '
    # remove mate information from read name and
    # concatenate CB and UMI with feature name
    cmd += ' awk \'{ sub(/\/[12]$/,"",$4); split($4,qname,/\//); $4=qname[2]"\\t"qname[3]"\\t"$16 } '
    #cmd += ' !x[$4]++ {OFS="\\t"; print $4}\' | '
    cmd += ' {OFS="\\t"; print $4}\' | '
    cmd += f' gzip > {isecFile}'

    writerr(f'Extracting {chrom} reference', verbose)
    run_shell_cmd(cmd0)
    writerr(f'Intersecting alignments with {chrom} reference', verbose)
    run_shell_cmd(cmd)
    writerr(f'Finished mapping {chrom}', verbose)

    return isecFile

# Concatenate and sort data obtained from isec()
def chrcat(filesList, threads, outdir, tmpdir, verbose):
    makedirs(outdir, exist_ok=True)
    mappings_file = tmpdir + '/cb_umi_te.bed.gz'
    barcodes_file = outdir + '/barcodes.tsv.gz'
    features_file = outdir + '/features.tsv.gz'
    bedFiles = ' '.join(filesList)
    cmd0 = f'zcat {bedFiles} '
    cmd0 += f' | LC_ALL=C sort --parallel {threads} --buffer-size 2G '
    cmd0 += f' | gzip > {mappings_file} '
    cmd1 = f'zcat {mappings_file} | cut -f1 | uniq | gzip > {barcodes_file} '
    cmd2 = f'zcat {mappings_file} '
    cmd2 += ' | awk \'!x[$3]++ { '
    #cmd2 += ' split($3,a,"~"); OFS="\\t"; print a[1],a[2],"Gene Expression" '
    cmd2 += ' split($3,a,"~"); '
    # avoid subfamilies with the same name
    cmd2 += ' if(a[1] in sf) { sf[a[1]]+=1 } else { sf[a[1]] }; '
    cmd2 += ' print a[1] sf[a[1]] "\\t" a[2] "\\tGene Expression" '
    cmd2 += ' }\' '
    cmd2 += f' | sort -u | gzip > {features_file} '

    writerr('Concatenating mappings', verbose)
    run_shell_cmd(cmd0)
    writerr(f'Writing mapped barcodes to {barcodes_file}')
    run_shell_cmd(cmd1)
    writerr(f'Writing mapped features to {features_file}')
    run_shell_cmd(cmd2)
    return mappings_file, barcodes_file, features_file
