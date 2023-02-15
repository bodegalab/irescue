# Dictionary with RepeatMasker files urls for each genome assembly, along with
# the number of lines taken by the header.

__genomes__ = {

    # human
    'GRCh38': ('https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.out.gz', 3),
    'hg38': ('https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.out.gz', 3),
    'GRCh37': ('https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.out.gz', 3),
    'hg19': ('https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.out.gz', 3),

    # mouse
    'GRCm39': ('https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.out.gz', 3),
    'mm39': ('https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.out.gz', 3),
    'GRCm38': ('https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/initial/mm10.fa.out.gz', 3),
    'mm10': ('https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/initial/mm10.fa.out.gz', 3),

    # fly
    'dm6': ('https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.out.gz', 3),

    # test
    'test': ('https://raw.githubusercontent.com/bodegalab/irescue/main/tests/data/rmsk.out.gz', 3)

}