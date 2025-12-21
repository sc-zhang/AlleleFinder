# AlleleFinder: a tool for identifying allele genes from polyploid genome.

[![DOI](https://zenodo.org/badge/284926541.svg)](https://doi.org/10.5281/zenodo.14015588)

## Introduction

This software is used for identifying allele genes from polyploid genome.

## Overview

![workflow](images/AlleleFinderWorkflow.png)

## Dependencies

Software:

- [MCScanX](https://github.com/wyp1125/MCScanX)
- [GMAP](http://research-pub.gene.com/gmap/)
- [NCBI BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)

> **Notice:** these software should be added to the <kbd>PATH</kbd> Environment Variable or create an environment use
> conda with environment.yaml.

## Installation

- Install latest release

```bash
# Download latest release
cd /path/to/install
wget https://github.com/sc-zhang/AlleleFinder/archive/refs/tags/V1.2.1.tar.gz
tar zxvf V1.2.1.tar.gz
chmod +x AlleleFinder-1.2.1/allelefinder.py
```

- Install development version

```bash
# Clone development repository
cd /path/to/install
git clone https://github.com/sc-zhang/AlleleFinder.git
chmod +x AlleleFinder/allelefinder.py
```

- Optional

```bash
# Optional, create environment with conda
conda env create -f environment.yaml
conda activate allelefinder_env

# Optional, add to PATH environment variable
echo 'export PATH=/path/to/install/AlleleFinder:$PATH' >> ~/.bash_profile
source ~/.bash_profile
```

## Usage

### 1. main program

```bash
usage: allelefinder.py [-h] [-v] {construct,cleanup,stat,adjust} ...

options:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit

Sub commands:
  {construct,cleanup,stat,adjust}
    construct           Construct allele table
    cleanup             Remove same CDS from cds, pep and gff3 files
    stat                Statistic allele table
    adjust              Adjust allele table with too many genes be marked as paralog
```

### 2. construct allele table

#### Usage

```bash
usage: allelefinder.py construct [-h] -r REF -d REF_CDS -f REF_GFF3 -c CDS -g GFF3 -n NUM_ALLELE [-m] [--ovlp_ratio OVLP_RATIO] [-b BLAST_ROUND] [--blast_threshold BLAST_THRESHOLD] [--reciprocal] [-e TE] [-j TE_OVERLAP] [--paralog_only] [-w WORKDIR]
                                 [-t THREADS]

options:
  -h, --help            show this help message and exit
  -r, --ref REF         reference fasta
  -d, --ref_cds REF_CDS
                        CDS fasta of ref
  -f, --ref_gff3 REF_GFF3
                        GFF3 file of ref
  -c, --cds CDS         CDS fasta of polyploid
  -g, --gff3 GFF3       GFF3 file of polyploid
  -n, --num_allele NUM_ALLELE
                        number of allele
  -m, --is_mono         if your reference fasta is mono assembly of polyploid, add this argument
  --ovlp_ratio OVLP_RATIO
                        threshold of gene pair coordinate overlap identified by GMAP, default: 0.8
  -b, --blast_round BLAST_ROUND
                        blast round, default: 2
  --blast_threshold BLAST_THRESHOLD
                        threshold of blast result which defined as matches*2/(query_length+reference_length), default: 0.5
  --reciprocal          add genes to allele table by blast with reciprocal, if set, blast_threshold would be ignored, default: False
  -e, --TE TE           TE gff3 for filtering, default: ""
  -j, --TE_overlap TE_OVERLAP
                        threshold of TE overlap, default: 0.3, only effect when TE is not NULL
  --paralog_only        do TE filter only on paralog genes
  -w, --workdir WORKDIR
                        workdir, default: wrkdir
  -t, --threads THREADS
                        threads, default: 12
```

> **Notice:**
> 1. the name of Chromosomes should be like: Chr01X, "X" means consecutive uppercase letters from A to Z, indicates
     different alleles, for example, if there are 4 alleles, the names should be: Chr01A,Chr01B,Chr01C,Chr01D.
> 2. the gff3 files must contain "gene" records, or you can use "sed" command to change "mRNA" to "gene" for some
     downloaded gff3 files.
> 3. there must no '-' in gene id.

#### Results

1. **Without TE filter**

   **allele.adjusted.txt** is the file contain all allele genes

   **allele.adjusted.*.stat** are the statistics information of allele

2. **With TE filter**

   **allele.adjusted.nonTEs.txt** is the file contain all allele genes

   **allele.adjusted.nonTEs.*.stat** are the statistics information of allele

### 3. Statistic

```bash
usage: allelefinder.py stat [-h] -i INPUT -g GFF3 -o OUTPUT

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input allele table
  -g GFF3, --gff3 GFF3  GFF3 file of polyploid
  -o OUTPUT, --output OUTPUT
                        Prefix of output file
```

### 4. Additional

* Cleanup sequence

If user want to remove the duplicate genes with same CDS/PEP sequence, the "cleanup" would deal with them and only one
gene
would be kept randomly.

```bash
usage: allelefinder.py cleanup [-h] [--in_cds IN_CDS] [--in_pep IN_PEP] --in_gff3 IN_GFF3 [--out_cds OUT_CDS] [--out_pep OUT_PEP] --out_gff3 OUT_GFF3 [--by_pep]

options:
  -h, --help           show this help message and exit
  --in_cds IN_CDS      Input CDS file
  --in_pep IN_PEP      Input PEP file
  --in_gff3 IN_GFF3    Input GFF3 file
  --out_cds OUT_CDS    Output CDS file
  --out_pep OUT_PEP    Output PEP file
  --out_gff3 OUT_GFF3  Output GFF3 file
  --by_pep             filter by PEP sequences
```

> **Notice:** CDS file and GFF3 file are required, PEP file is optional.

* Rescue genes which be cleaned up

If the cleanup method was applied, user may want to add genes which be cleanup onto allele table again the rescue method
would to do like this.

```bash
usage: allelefinder.py rescue [-h] -i INPUT --gff3 GFF3 --cds CDS [--pep PEP] -o OUTPUT

options:
  -h, --help           show this help message and exit
  -i, --input INPUT    Input allele table
  --gff3 GFF3          Input GFF3 file before cleanup
  --cds CDS            Input CDS file before cleanup
  --pep PEP            Input PEP file before cleanup, required when cleanup was run with --by_pep
  -o, --output OUTPUT  Output rescued allele table
```

* Adjust paralog genes

If there are too many genes be marked with paralog, you can use command below to pull them down as new alleles

```bash
usage: allelefinder.py adjust [-h] -i INPUT -m MIN_NUM -o OUTPUT

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input allele table
  -m MIN_NUM, --min_num MIN_NUM
                        Minium number of genes, which means the number of genes marked as paralog that distribute in different allele should be pulled down as new allele genes
  -o OUTPUT, --output OUTPUT
                        Output allele table
```
