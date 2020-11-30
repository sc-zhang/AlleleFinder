#!/usr/bin/env python
import sys
import re


def adjust_allele_table(in_allele, mono_gff3, hap_gff3, blast, iden, tandem, allele_num, out_allele):
	print("Loading blast")
	blast_db = {}
	with open(blast, 'r') as fin:
		for line in fin:
			data = line.strip().split()
			hap_gn = data[0]
			mono_gn = data[1]
			bi = float(data[2])
			if bi < iden:
				continue
			bs = data[-1]
			if hap_gn not in blast_db:
				blast_db[hap_gn] = [mono_gn, bs]
			if bs > blast_db[hap_gn][-1]:
				blast_db[hap_gn] = [mono_gn, bs]
	
	print("Loading mono gff3")
	mono_db = {}
	with open(mono_gff3, 'r') as fin:
		for line in fin:
			if line.strip() == '' or line[0] == '#':
				continue
			data = line.strip().split()
			if data[2] == 'gene':
				id = data[8].split(";")[0].split('=')[1]
				if id not in mono_db:
					mono_db[id] = [data[0], int(data[3])]
	
	print("Loading hap_gff3")
	hap_db = {}
	with open(hap_gff3, 'r') as fin:
		for line in fin:
			if line.strip() == '' or line[0] == '#':
				continue
			data = line.strip().split()
			if data[2] != 'gene':
				continue
			chrn = data[0]
			id = data[8].split(";")[0].split('=')[1]
			hap_db[id] = [chrn, int(data[3])]
	
	print("Loading tandem")
	tandem_db = {}
	with open(tandem, 'r') as fin:
		for line in fin:
			data = line.strip().split(',')
			for id in data:
				tandem_db[id] = 1
	
	print("Formatting allele table")
	full_allele = []
	info_db = {}
	with open(in_allele, 'r') as fin:
		for line in fin:
			data = line.strip().split(',')
			tmp_db = {}
			for id in data:
				hchrn, hpos = hap_db[id]
				hidx = ord(hchrn[-1])-65
				if id not in blast_db:
					continue
				mono_gn = blast_db[id][0]
				if mono_gn not in mono_db:
					continue
				if mono_gn not in info_db:
					info_db[mono_gn] = [[] for i in range(0, allele_num)]
				info_db[mono_gn][hidx].append([hpos, id])
		
		for mono_gn in info_db:
			tmp_list = []
			null_cnt = 0
			for i in range(0, len(info_db[mono_gn])):
				if info_db[mono_gn][i] == []:
					info_db[mono_gn][i] = "NA"
					null_cnt += 1
				else:
					tmp_list = sorted(info_db[mono_gn][i])
					tmp_ids = []
					tmp_ids.append(tmp_list[0][-1])
					if len(tmp_list) > 1:
						for j in range(1, len(tmp_list)):
							pos, id = tmp_list[j]
							if id in tandem_db:
								tmp_ids.append(id+"-T")
							else:
								tmp_ids.append(id+"-P")
					info_db[mono_gn][i] = ','.join(tmp_ids)

			if null_cnt == allele_num:
				continue
			
			tmp_list = mono_db[mono_gn]
			tmp_list.append(mono_gn)
			tmp_list.extend(info_db[mono_gn])
			full_allele.append(tmp_list)
			
	
	print("Writing allele table")
	with open(out_allele, 'w') as fout:
		with open(out_allele+'.simple', 'w') as fsimple:
			allele_header = []
			for i in range(0, allele_num):
				allele_header.append("Allele %s"%(chr(i+65)))
			fout.write("#CHR\tPOS\tMono gene\t%s\n"%('\t'.join(allele_header)))
			cnt = 0
			for info in sorted(full_allele):
				fout.write("%s\n"%('\t'.join(map(str, info))))
				tmp_list = []
				for id in info[3:]:
					if id != "NA":
						tmp_list.append(id)
				fsimple.write("%s\n"%(','.join(tmp_list).replace('-T', '').replace('-P', '')))
	print("Finished")


if __name__ == "__main__":
	if len(sys.argv) < 9:
		print("Usage: python %s <in_allele> <mono_gff3> <hap_gff3> <hap_mono_blast> <iden_threshold> <tandem_list> <allele_number> <out_allele>"%sys.argv[0])
	else:
		in_allele, mono_gff3, hap_gff3, blast, iden, tandem, allele_num, out_allele = sys.argv[1:]
		adjust_allele_table(in_allele, mono_gff3, hap_gff3, blast, float(iden), tandem, int(allele_num), out_allele)
