![GitHub Workflow Status](https://img.shields.io/github/actions/workflow/status/bodegalab/irescue/python-publish.yml?logo=github&label=build)
[![PyPI](https://img.shields.io/pypi/v/irescue?logo=python)](https://pypi.org/project/irescue/)
[![container](https://img.shields.io/badge/dynamic/json?url=https%3A%2F%2Fquay.io%2Fapi%2Fv1%2Frepository%2Fbiocontainers%2Firescue%2Ftag%2F&query=%24.tags.0.name&logo=docker&label=docker%2Fsingularity&color=%231D63ED)](#container)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat&logo=anaconda)](https://bioconda.github.io/recipes/irescue/README.html)
[![paper](https://img.shields.io/badge/Nucleic%20Acids%20Res-10.1093%2Fnar%2Fgkae793-orange)](https://doi.org/10.1093/nar/gkae793)
[![zenodo](https://img.shields.io/badge/Zenodo-10.5281/zenodo.13479363-blue)](https://doi.org/10.5281/zenodo.13479363)

# IRescue - <ins>I</ins>nterspersed <ins>Re</ins>peats <ins>s</ins>ingle-<ins>c</ins>ell q<ins>u</ins>antifi<ins>e</ins>r

<img align="right" height="160" src="docs/logo.png">
IRescue quantifies the expression fo transposable elements (TEs) subfamilies in single cell RNA sequencing (scRNA-seq) data, performing UMI-deduplication with sequencing errors correction and probabilistic assignment of multi-mapping reads by Expectation-Maximization (EM). The output is written on a sparse matrix compatible with Seurat, Scanpy and other toolkits.

## Content

- [Installation](#installation)
  - [Using conda](#conda)
  - [Using pip](#pip)
  - [Container (Docker/Singularity)](#container)
- [Usage](#usage)
  - [Quick start](#quick_start)
  - [Output files](#output_files)
  - [Load IRescue data with Seurat](#seurat)
- [Cite](#cite)

## <a name="installation"></a>Installation

### <a name="conda"></a>Using conda (recommended)

We recommend using conda, as it will install all the required packages along IRescue.

```bash
conda create -n irescue -c conda-forge -c bioconda irescue
```

### <a name="pip"></a>Using pip

If for any reason it is not possible or desired to use conda, it can be installed with pip and the following requirements must be installed manually: `python>=3.9`, `samtools>=1.12`, `bedtools>=2.30.0`, and fairly recent versions of the GNU utilities are required, e.g. `gawk>=5.0.1`, `coreutils>=8.30` and `gzip>=1.10` (older versions are untested).

```bash
pip install irescue
```

### Build from source

By building the package directly from the source, you can try out the features and bug fixes that will be implemented in the future release. As above, you need to install some requirements manually. Be aware that builds from the development branches may be unstable.

```bash
git clone https://github.com/bodegalab/irescue
cd irescue
pip install .
```

### <a name="container"></a>Container (Docker/Singularity)

Docker and Singularity containers are available for each conda release of IRescue. Choose the `TAG` corresponding to the desired IRescue version [from the Biocontainers repository](https://quay.io/repository/biocontainers/irescue?tab=tags) and pull or execute the container with Docker or Singularity:

```bash
# Get latest biocontainers tag (with curl and python3, otherwise check the above link for the desired version/tag)
TAG=$(curl -s -X GET https://quay.io/api/v1/repository/biocontainers/irescue/tag/ | python3 -c 'import json,sys;obj=json.load(sys.stdin);print(obj["tags"][0]["name"])')

# Run with Docker
docker run quay.io/biocontainers/irescue:$TAG irescue --help

# Run with Singularity
singularity exec https://depot.galaxyproject.org/singularity/irescue:$TAG irescue --help
```

## <a name="usage"></a>Usage

```sh
irescue --help
```

The only required input is a BAM file annotated with cell barcode and UMI sequences as tags (by default, `CB` tag for cell barcode and `UR` tag for UMI; override with `--cb-tag` and `--umi-tag`).

You can obtain it by aligning your reads using [STARsolo](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md). It is advised to keep secondary alignments in BAM file, that will be used in the EM procedure to assign multi-mapping reads (e.g. `--outFilterMultimapNmax 100 --winAnchorMultimapNmax 100` or more), and remember to output all the needed SAM attributes (e.g. `--outSAMattributes NH HI AS nM NM MD jM jI XS MC ch cN CR CY UR UY GX GN CB UB sM sS sQ`).

RepeatMasker annotation will be automatically downloaded for the chosen genome assembly (e.g. `-g hg38`), or provide your own annotation in bed format (e.g. `-r TE.bed`).

```sh
irescue -b genome_alignments.bam -g hg38
```

If you already obtained gene-level counts (using STARsolo, Cell Ranger, Alevin, Kallisto or other tools), it is advised to provide the whitelisted cell barcodes list as a text file (`-w barcodes.tsv`). This will significantly improve performance by processing viable cells only.

For optimal run time, use at least, e.g.: `-p 8`.

### <a name="output_files"></a>Output files

IRescue generates TE counts in a sparse matrix readable by [Seurat](https://github.com/satijalab/seurat) or [Scanpy](https://github.com/scverse/scanpy) into a `counts/` subdirectory. Optional outputs include a description of equivalence classes with UMI deduplication stats `ec_dump.tsv.gz` and a subdirectory of temporary files `tmp/` for debugging purpose. A detailed logging is enabled by `--verbose` and written to standard error.

```
irescue_out/
├── counts/
│   ├── barcodes.tsv.gz
│   ├── features.tsv.gz
│   └── matrix.mtx.gz
├── ec_dump.tsv.gz
└── tmp/
```

### <a name="seurat"></a>Load IRescue data with Seurat

To integrate TE counts into an existing Seurat object containing gene expression data, they can be added as an additional assay:

```R
# import TE counts from IRescue output directory
te.data <- Seurat::Read10X('./IRescue_out/', gene.column = 1, cell.column = 1)

# create Seurat assay from TE counts
te.assay <- Seurat::CreateAssayObject(te.data)

# subset the assay by the cells already present in the Seurat object (in case it has been filtered)
te.assay <- subset(te.assay, colnames(te.assay)[which(colnames(te.assay) %in% colnames(seurat_object))])

# add the assay in the Seurat object
seurat_object[['TE']] <- irescue.assay
```

The result will be something like this:
```
An object of class Seurat 
32276 features across 42513 samples within 2 assays 
Active assay: RNA (31078 features, 0 variable features)
 1 other assay present: TE
```

From here, TE expression can be normalized. To normalize according to gene counts or TE+gene counts, normalize manually or merge the assays. Reductions can be made using TE, gene or TE+gene expression.

## <a name="cite"></a>Cite

Benedetto Polimeni, Federica Marasca, Valeria Ranzani, Beatrice Bodega, **IRescue: uncertainty-aware quantification of transposable elements expression at single cell level**, *Nucleic Acids Research*, 2024; https://doi.org/10.1093/nar/gkae793
