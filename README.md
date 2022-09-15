# IRescue

IRescue is a software for quantifying the expression of transposable elements (TEs) subfamilies in single cell RNA sequencing (scRNA-seq) data.
The core feature of IRescue is to consider all multiple alignments (i.e. non-primary alignments) of reads/UMIs mapping on multiple TEs in a BAM file, to better infer the TE subfamily of origin. IRescue implements a UMI error-correction, deduplication and quantification strategy that includes such alignment events. IRescue's output is compatible with most scRNA-seq analysis toolkits, such as Seurat or Scanpy.

## Content

- [Installation](#installation)
- [Usage](#usage)
  - [Quick start](#quick_start)
  - [Output files](#output_files)
  - [Load IRescue data with Seurat](#seurat)

## <a name="installation"></a>Installation

We advise installing IRescue in a conda environment with all required software:

```bash
# clone this repository
git clone https://github.com/bepoli/irescue
cd irescue

# create a conda environment
conda env create -f docker/environment.yaml -n irescue
conda activate irescue

# install irescue
pip install .

# test the installation
irescue --help
```

## <a name="usage"></a>Usage

### <a name="quick_start"></a>Quick start

The only required input is a BAM file annotated with cell barcode and UMI sequences as tags (by default, `CB` tag for cell barcode and `UR` tag for UMI; override with `--CBtag` and `--UMItag`). You can obtain it by aligning your reads using [STARsolo](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md).

RepeatMasker annotation will be automatically downloaded for the chosen genome assembly (e.g. `-g hg38`), or provide your own annotation in bed format (e.g. `-r TE.bed`).

```bash
irescue -b genome_alignments.bam -g hg38
```

If you already obtained gene-level counts (using STARsolo, Cell Ranger, Alevin, Kallisto or other tools), it is advised to provide the whitelisted cell barcodes list as a text file, e.g.: `-w barcodes.tsv`. This will significantly improve performance.

IRescue performs best using at least 4 threads, e.g.: `-p 8`.

### <a name="output_files"></a>Output files

IRescue generates TE counts in a sparse matrix format, readable by [Seurat](https://github.com/satijalab/seurat) or [Scanpy](https://github.com/scverse/scanpy):

```
IRescue_out/
├── barcodes.tsv.gz
├── features.tsv.gz
└── matrix.mtx.gz
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
