def allele_blast(in_blast, iden):
	allele_list = []
	used_genes = {}
	with open(in_blast, 'r') as fin:
		for line in fin:
			data = line.strip().split()
			qry = data[0]
			target = data[1]
			if float(data[2]) < iden:
				continue
			if qry not in used_genes:
				used_genes[qry] = 1
				allele_list.append([qry, target])
	return allele_list
