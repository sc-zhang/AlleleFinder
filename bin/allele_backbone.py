import re


def merge_allele(allele_list):
	pos_db = {}
	new_allele_list = []
	for i in range(0, len(allele_list)):
		merge_list = []
		scrappy_allele = []
		for allele in allele_list[i]:
			if allele in pos_db:
				merge_list.append(allele)
			else:
				scrappy_allele.append(allele)
		if merge_list == []:
			new_allele_list.append(allele_list[i])
			for allele in allele_list[i]:
				pos_db[allele] = len(new_allele_list)-1
		else:
			merge_pos = pos_db[merge_list[0]]
			if scrappy_allele != []:
				for allele in scrappy_allele:
					pos_db[allele] = merge_pos
				new_allele_list[merge_pos].extend(scrappy_allele)
			for allele in merge_list[1:]:
				allele_pos = pos_db[allele]
				if allele_pos == merge_pos:
					continue
				new_allele_list[merge_pos].extend(new_allele_list[allele_pos])
				for other_allele in new_allele_list[allele_pos]:
					pos_db[other_allele] = merge_pos
				new_allele_list[allele_pos] = []
	allele_list = []
	for i in range(0, len(new_allele_list)):
		if new_allele_list[i] != []:
			allele_list.append(list(set(new_allele_list[i])))
	
	return allele_list


def split_allele(allele_list, gene_db):
	new_allele_list = []
	for i in range(0, len(allele_list)):
		tmp_list = allele_list[i]
		tmp_db = {}
		for gene in tmp_list:
			if gene not in gene_db:
				continue
			chrn = gene_db[gene][0]
			if chrn not in tmp_db:
				tmp_db[chrn] = []
			tmp_db[chrn].append(gene)
		for chrn in tmp_db:
			new_allele_list.append(tmp_db[chrn])			
	return new_allele_list


def get_allele_with_mcscanx(html_file):
	allele_list = []
	with open(html_file, 'r') as fin:
		for line in fin:
			if line.startswith("<tr"):
				vals = re.findall('<td.*?>(.*?)</td>', line)
				tmp_list = []
				if len(vals) < 2:
					continue
				for val in vals[1:]:
					if '&' not in val and 'ctg' not in val:
						tmp_list.append(val)
				if tmp_list != []:
					allele_list.append(tmp_list)
	return allele_list
