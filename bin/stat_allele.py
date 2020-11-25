import re


def stat_allele(in_csv, in_gff3, allele_num, out_stat):
	gene_db = {}
	chrn_list = {}
	with open(in_gff3, 'r') as fin:
		for line in fin:
			if line.strip() == '' or line[0] == '#':
				continue
			data = line.strip().split()
			if "Name" in data[8]:
				id = re.findall(r'Name=(.*)', data[8])[0].split(';')[0]
			else:
				id = re.findall(r'ID=(.*)', data[8])[0].split(';')[0]
			if 'tig' in data[0] or 'ctg' in data[0]:
				continue
			chridx = re.findall(r'[A-Za-z]*(\d+)', data[0])[0]
			chrn = "Chr%02d"%(int(chridx))
			if chrn not in chrn_list:
				chrn_list[chrn] = data[0]
			gene_db[id] = chrn
	
	allele_db = {}
	gene_cnt = {}
	with open(in_csv, 'r') as fin:
		for line in fin:
			data = line.strip().split(',')
			chrn = gene_db[data[0]]
			if chrn not in allele_db:
				gene_cnt[chrn] = 0
				allele_db[chrn] = {}
				for i in range(1, allele_num+2):
					allele_db[chrn][i] = 0
			allele_cnt = 1
			gene_cnt[chrn] += 1
			for gene in data[1:]:
				if gene.strip() != '':
					allele_cnt += 1
			if allele_cnt <= allele_num:
				allele_db[chrn][allele_cnt] += 1
			else:
				allele_db[chrn][allele_num+1] += 1
	with open(out_stat, 'w') as fout:
		out_info = ['', 'Total', 'No. more than %d'%allele_num]
		total_db = {}
		total_gene_cnt = 0
		for i in range(allele_num, 0, -1):
			out_info.append("No. with %d"%i)
		fout.write(','.join(out_info)+'\n')
		for chrn in sorted(allele_db):
			out_info = [chrn, str(gene_cnt[chrn])]
			for cnt in sorted(allele_db[chrn], reverse=True):
				out_info.append(str(allele_db[chrn][cnt]))
				if cnt not in total_db:
					total_db[cnt] = 0
				total_db[cnt] += allele_db[chrn][cnt]
			total_gene_cnt += gene_cnt[chrn]
			fout.write(','.join(out_info)+'\n')
		out_info = ['Total', str(total_gene_cnt)]
		for cnt in sorted(total_db, reverse=True):
			out_info.append(str(total_db[cnt]))
		fout.write(','.join(out_info)+'\n')