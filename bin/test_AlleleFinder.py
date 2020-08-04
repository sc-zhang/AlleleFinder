#!/usr/bin/env python
import sys
import allele_backbone as ab
import allele_gmap as ag
import allele_blast as abl

in_html = "/public6/zsc/AlleleFinder/test/xyz/xyz.html/Ss.html"
in_gff3 = "/public6/zsc/AlleleFinder/test/gmap.gff"
in_blast = "/public6/zsc/AlleleFinder/test/single.blast"
in_blast2 = "/public6/zsc/AlleleFinder/test/single2.blast"


base_allele = ab.get_allele_with_mcscanx(in_html)
gff3_db, gene_order = ag.read_gff3(in_gff3)
gff3_allele = ag.allele_gmap(gff3_db, 12)
blast_allele = abl.allele_blast(in_blast)
blast2_allele = abl.allele_blast(in_blast2)

base_allele.extend(gff3_allele)

base_allele = ab.split_allele(base_allele)

base_allele.extend(blast_allele)
base_allele.extend(blast2_allele)

allele_list = ab.merge_allele(base_allele)

with open("origin.csv", 'w') as fout:
	for allele in sorted(allele_list):
		fout.write(",".join(sorted(allele)))
		fout.write("\n")
