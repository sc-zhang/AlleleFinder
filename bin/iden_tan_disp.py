import re


def read_blast(in_blast):
	blast_db = {}
	with open(in_blast, 'r') as fin:
		for line in fin:
			data = line.strip().split()
			qry = data[0]
			target = data[1]
			if qry not in blast_db:
				blast_db[qry] = []
			blast_db[qry].append(target)
	return blast_db


def get_tandam(allele_list, gene_order, blast_db):
	for allele in allele_list:
		inter = set(blast_db[allele[0]])
		for gene in allele:
			inter = inter.intersection(set(blast_db[gene]))
		match = ""
		if len(inter) != 0:
			for target in blast_db[gene]:
				if target in inter:
					match = target
					break
		chr_cnt = {}
		for gene in allele:
			pre, chridx, ga, gn = re.findall(r'(\S+)\.(\d+)([A-Za-z])(\d+)', gene)
			if chridx not in chr_cnt:
				chr_cnt[chridx] = 0
			chr_cnt[chridx] += 1
