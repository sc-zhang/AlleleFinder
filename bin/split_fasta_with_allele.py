#!/usr/bin/env python
import sys


def split_fasta_with_allele(in_fa, in_allele, out_pre):
	fa_db = {}
	with open(in_fa, 'r') as fin:
		for line in fin:
			if line[0] == '>':
				id = line.strip().split()[0][1:]
				fa_db[id] = []
			else:
				fa_db[id].append(line.strip())
	
	for id in fa_db:
		fa_db[id] = ''.join(fa_db[id])
	
	single_list = []
	multi_list = []
	with open(in_allele, 'r') as fin:
		for line in fin:
			data = line.strip().split(',')
			cnt = 0
			for gene in data:
				if gene.strip() != '':
					cnt += 1
			if cnt == 1:
				single_list.append(data[0])
			else:
				for gene in data:
					if gene.strip() != '':
						multi_list.append(gene)
		
	with open(out_pre+'_single.fa', 'w') as fout:
		for id in single_list:
			fout.write(">%s\n%s\n"%(id, fa_db[id]))
	
	with open(out_pre+'_multi.fa', 'w') as fout:
		for id in multi_list:
			fout.write(">%s\n%s\n"%(id, fa_db[id]))


if __name__ == "__main__":
	if len(sys.argv) < 4:
		print("Usage: python %s <in_fasta> <in_allele> <out_pre>"%sys.argv[0])
	else:
		in_fa, in_allele, out_pre = sys.argv[1:]
		split_fasta_with_allele(in_fa, in_allele, out_pre)

