## Introduction

This software is used for identifying allele genes from polyploid genome.



## Dependencies

Software:

&ensp;&ensp;&ensp;&ensp;MCScanX  
&ensp;&ensp;&ensp;&ensp;GMAP  
&ensp;&ensp;&ensp;&ensp;NCBI BLAST+  



## Installation

```bash
cd /path/to/install
git clone https://github.com/sc-zhang/AlleleFinder.git
chmod +x AlleleFinder/bin/AlleleFinder
# Optional
echo 'export PATH=/path/to/install/AlleleFinder/bin:$PATH' >> ~/.bash_profile
source ~/.bash_profile
```



## Usage

```bash
usage: AlleleFinder [-h] -r REF -d REF_CDS -f REF_GFF3 -c CDS -g GFF3 -n NUM_ALLELE [-m] [-b BLAST_COUNT] [-i BLAST_IDENTITY] [-e TE] [-j TE_OVERLAP] [-w WORKDIR] [-t THREADS]

options:
  -h, --help            show this help message and exit
  -r REF, --ref REF     reference fasta
  -d REF_CDS, --ref_cds REF_CDS
                        CDS fasta of ref
  -f REF_GFF3, --ref_gff3 REF_GFF3
                        GFF3 file of ref
  -c CDS, --cds CDS     CDS fasta of polyploid
  -g GFF3, --gff3 GFF3  GFF3 file of polyploid
  -n NUM_ALLELE, --num_allele NUM_ALLELE
                        number of allele
  -m, --is_mono         If your reference fasta is mono assembly of polyploid, add this argument
  -b BLAST_COUNT, --blast_count BLAST_COUNT
                        blast count, default: 2
  -i BLAST_IDENTITY, --blast_identity BLAST_IDENTITY
                        threshold of blast identity, default: 80
  -e TE, --TE TE        TE gff3 for filtering, default: ""
  -j TE_OVERLAP, --TE_overlap TE_OVERLAP
                        threshold of TE overlap, default: 0.3, only effect when TE is not NULL
  -w WORKDIR, --workdir WORKDIR
                        workdir, default: wrkdir
  -t THREADS, --threads THREADS
                        threads, default: 12
```

**Notice:** the name of Chromosomes should be like: Chr01X, "X" means consecutive uppercase letters from A to Z, means 
different alleles, for example, if there are 4 alleles, the names should be: Chr01A,Chr01B,Chr01C,Chr01D  
**Notice:** the gff3 files must contain "gene" records, or you can use "sed" command to change "mRNA" to "gene" for some downloaded gff3 files.  
**Notice:** there must no '-' in gene id


## Results

1. **Without TE filter**

**allele.adjusted.txt** is the file contain all allele genes

**allele.adjusted.*.stat** are the statistics information of allele

2. **With TE filter**

**allele.adjusted.nonTEs.txt** is the file contain all allele genes

**allele.adjusted.nonTEs.*.stat** are the statistics information of allele


## Additional
If there are too many genes be marked with paralog, you can use pull_down_paralog.py to pull them down as new alleles
```bash
usage: pull_down_paralog.py [-h] -i INPUT -m MIN_NUM -o OUTPUT

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input allele table
  -m MIN_NUM, --min_num MIN_NUM
                        Minium number of genes, which means the number of genes marked as paralog that distribute in different allele should be pulled down as new allele genes
  -o OUTPUT, --output OUTPUT
                        Output allele table
```
**Notice** because we only pull down the first paralog genes from each allele to contruct new allele, that means if there are more than one paralog genes in different alleles, you may need run this script more than one time to pull all paralog genes which with the distribution mentioned before down as new alleles.