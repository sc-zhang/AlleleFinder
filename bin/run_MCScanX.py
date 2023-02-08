#!/usr/bin/env python
import sys
import os
import re


def check_fasta_type(in_fa):
	with open(in_fa, 'r') as fin:
		for line in fin:
			if line[0] == ">":
				continue
			if "M" in line.upper():
				return "PROT"
			else:
				return "NUCL"


def run_MCScanX(in_fa, in_gff3, out_dir, threads):
	in_fa = os.path.abspath(in_fa)
	in_gff3 = os.path.abspath(in_gff3)
	if not os.path.exists(out_dir):
		os.mkdir(out_dir)
	
	print("Entering %s"%out_dir)
	os.chdir(out_dir)

	if not os.path.exists("xyz"):
		os.mkdir("xyz")
	
	if not os.path.exists("xyz/xyz.blast"):
		print("Running blast")
		fa_type = check_fasta_type(in_fa)
		if fa_type == "NUCL":
			cmd1 = "makeblastdb -in %s -dbtype nucl -out blastdb &> /dev/null"%in_fa
			cmd2 = "blastn -query %s -db blastdb -out xyz/xyz.blast -evalue 1e-10 -outfmt 6 -num_alignments 5 -num_threads %d &> /dev/null"%(in_fa, threads)
		else:
			cmd1 = "makeblastdb -in %s -dbtype prot -out blastdb &> /dev/null"%in_fa
			cmd2 = "blastp -query %s -db blastdb -out xyz/xyz.blast -evalue 1e-10 -outfmt 6 -num_alignments 5 -num_threads %d &> /dev/null"%(in_fa, threads)
		
		print("Running command: %s"%cmd1)
		os.system(cmd1)
		print("Running command: %s"%cmd2)
		os.system(cmd2)
	else:
		print("Blast file found, skip")

	if not os.path.exists("xyz/xyz.gff"):
		print("Loading gff3")
		gff3_db = {}
		with open(in_gff3, 'r') as fin:
			for line in fin:
				if line.strip() == '' or line[0] == '#':
					continue
				data = line.strip().split()
				chrn = re.findall(r'([A-Za-z]+\d+)', data[0])[0]
				if 'tig' in chrn or 'ctg' in chrn:
					continue
				if data[2] == 'gene':
					if "Name" in data[8]:
						id = re.findall(r'Name=(.*)', data[8])[0].split(';')[0]
					else:
						id = re.findall(r'ID=(.*)', data[8])[0].split(';')[0]
					
					if chrn not in gff3_db:
						gff3_db[chrn] = []
					gff3_db[chrn].append([id, int(data[3]), int(data[4])])
		
		print("Writing gff")
		idx = 1
		with open(os.path.join("xyz/xyz.gff"), 'w') as fout:
			for chrn in sorted(gff3_db):
				spid = "NN%02d"%(idx)
				idx += 1
				for id, sp, ep in sorted(gff3_db[chrn]):
					fout.write("%s\t%s\t%d\t%d\n"%(spid, id, sp, ep))
	else:
		print("GFF file found, skip")
	
	cmd = "MCScanX xyz/xyz &> /dev/null"
	cmd_class = "duplicate_gene_classifier xyz/xyz"
	print("Running command: %s"%cmd)
	os.system(cmd)
	print("Running command: %s"%cmd_class)
	os.system(cmd_class)
	print("Finished")	


if __name__ == "__main__":
	if len(sys.argv) < 5:
		print("Usage: python %s <in_fasta> <in_gff3> <out_dir> <threads>"%sys.argv[0])
	else:
		in_fa, in_gff3, out_dir, threads = sys.argv[1:]
		run_MCScanX(in_fa, in_gff3, out_dir, int(threads))
