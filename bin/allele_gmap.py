import multiprocessing
import re


def read_gff3(in_gff3):
	gff3_db = {}
	gene_order = {}
	gene_cnt = {}
	with open(in_gff3, 'r') as fin:
		for line in fin:
			if line.strip() == '' or line[0] == '#':
				continue
			data = line.strip().split()
			chrn = data[0]
			if 'ctg' in chrn or 'tig' in chrn:
				continue
			type = data[2]
			if type != 'gene':
				continue
			sp = int(data[3])
			ep = int(data[4])
			if 'Name' in data[8]:
				gene = re.findall(r'Name=(.*)', data[8])[0].split(';')[0]
			else:
				gene = re.findall(r'ID=(.*)', data[8])[0].split(';')[0]
			if 'tig' in gene or 'ctg' in gene:
				continue
			if chrn not in gff3_db:
				gff3_db[chrn] = []
			if chrn not in gene_cnt:
				gene_cnt[chrn] = 0
			gene_cnt[chrn] += 1
			gene_order[gene] = [chrn, gene_cnt[chrn]]
			gff3_db[chrn].append([sp, ep, gene])
	return gff3_db, gene_order


def merge_allele(region_list):
	alleles = []
	last_ep = 0
	for sp, ep, gene in sorted(region_list):
		if last_ep == 0:
			alleles.append([])
			alleles[-1].append(gene)
			last_ep = ep
		else:
			if sp > last_ep:
				alleles.append([])
			alleles[-1].append(gene)
			last_ep = ep
	return alleles


def allele_gmap(gff3_db, threads):
	pool = multiprocessing.Pool(processes=threads)
	res = {}
	for chrn in gff3_db:
		res[chrn] = pool.apply_async(merge_allele, (gff3_db[chrn],))
	pool.close()
	pool.join()
	allele_list = []
	for chrn in res:
		allele_list.extend(res[chrn].get())
	
	return allele_list
