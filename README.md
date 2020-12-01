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
AlleleFinder [-h] -m MONO -d MONO_CDS -f MONO_GFF3 -c CDS -g GFF3 -n NUM_ALLELE [-b BLAST_COUNT] [-i BLAST_IDENTITY] [-e TE] [-j TE_OVERLAP] [-w WORKDIR] [-t THREADS]
```

**Notice** there must no '-' in gene id

**-m, --mono** fasta file of mono genome

**-d, --mono_cds** cds file of mono genome

**-f, --mono_gff3** gff3 file of mono genome

**-c, --cds** cds file of polyploid genome

**-g, --gff3** gff3 file of polyploid  genome

**-n, --num_allele** number of allele

**-b, --blast_count** iteration count for running blast, default: 2

**-i, --blast_identity** threshold of blast identity, default: 80

**-e, --TE** TE file for filtering, default: NULL

**-j, --TE_overlap** threshold of TE overlap, default: 0.3, only effect when TE is not NULL

**-w, --workdir** work directory, default: wrkdir

**-t, --threads** threads, default: 12



## Results

1. **Without TE filter**

**allele.adjusted.txt** is the file contain all allele genes

**allele.adjusted.*.stat** are the statistics information of allele

2. **With TE filter**

**allele.adjusted.nonTEs.txt** is the file contain all allele genes

**allele.adjusted.nonTEs.*.stat** are the statistics information of allele
