#!/usr/bin/env python
import sys
import re
import copy


def adjust_allele_table(in_allele, ref_gff3, hap_gff3, ref_blast, hap_blast, iden, tandem, allele_num, out_allele, is_mono):
	print("Loading ref blast")
	ref_blast_db = {}
	with open(ref_blast, 'r') as fin:
		for line in fin:
			data = line.strip().split()
			hap_gn = data[0]
			ref_gn = data[1]
			bi = float(data[2])
			if bi < iden:
				continue
			bs = data[-1]
			if hap_gn not in ref_blast_db:
				ref_blast_db[hap_gn] = [ref_gn, bs]
			if bs > ref_blast_db[hap_gn][-1]:
				ref_blast_db[hap_gn] = [ref_gn, bs]
	
	print("Loading hap blast")
	hap_blast_db = {}
	with open(hap_blast, 'r') as fin:
		for line in fin:
			data = line.strip().split()
			gn1 = data[0]
			gn2 = data[1]
			bi = float(data[2])
			if gn1 != gn2 and bi >= iden:
				if gn1 not in hap_blast_db:
					hap_blast_db[gn1] = [gn2, bi]
				elif bi > hap_blast_db[gn1][1]:
					hap_blast_db[gn1] = [gn2, bi]

	print("Loading ref gff3")
	ref_db = {}
	with open(ref_gff3, 'r') as fin:
		for line in fin:
			if line.strip() == '' or line[0] == '#':
				continue
			data = line.strip().split()
			if data[2] == 'gene':
				if "Name" in data[8]:
					id = re.findall(r'Name=(.*)', data[8])[0].split(';')[0]
				else:
					id = re.findall(r'ID=(.*)', data[8])[0].split(';')[0]
				if id not in ref_db:
					ref_db[id] = [data[0], int(data[3])]
	
	print("Loading hap gff3")
	hap_db = {}
	with open(hap_gff3, 'r') as fin:
		for line in fin:
			if line.strip() == '' or line[0] == '#':
				continue
			data = line.strip().split()
			if data[2] != 'gene':
				continue
			chrn = data[0]
			if "Name" in data[8]:
				id = re.findall(r'Name=(.*)', data[8])[0].split(';')[0]
			else:
				id = re.findall(r'ID=(.*)', data[8])[0].split(';')[0]
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
	ref_cnt = {}
	with open(in_allele, 'r') as fin:
		for line in fin:
			data = line.strip().split(',')
			ref_gn_cnt = {}
			tmp_list = [[] for i in range(0, allele_num)]
			for id in data:
				hchrn, hpos = hap_db[id]
				hidx = ord(hchrn[-1])-65
				if id not in ref_blast_db:
					if id not in hap_blast_db:
						continue
					ref_id = hap_blast_db[id][0]
				else:
					ref_id = id
				if ref_id not in ref_blast_db or ref_id not in data:
					continue
				ref_gn = ref_blast_db[ref_id][0]
				if ref_gn not in ref_db:
					continue
				if ref_gn not in ref_gn_cnt:
					ref_gn_cnt[ref_gn] = 0
				ref_gn_cnt[ref_gn] += 1
				if is_mono.upper() == 'T':
					if ref_gn not in info_db:
						info_db[ref_gn] = [[] for i in range(0, allele_num)]
					info_db[ref_gn][hidx].append([hpos, id])
				else:
					tmp_list[hidx].append([hpos, id])
			if is_mono.upper() == 'F':
				max_gn = ""
				max_cnt = 0
				for ref_gn in ref_gn_cnt:
					if ref_gn_cnt[ref_gn] > max_cnt:
						max_cnt = ref_gn_cnt[ref_gn]
						max_gn = ref_gn
				if max_gn not in ref_cnt:
					ref_cnt[max_gn] = 0
				ref_cnt[max_gn] += 1
				new_gn = "%s:::%d"%(max_gn, ref_cnt[max_gn])
				info_db[new_gn] = tmp_list
		
		for ref_gn in info_db:
			tmp_list = []
			null_cnt = 0
			for i in range(0, len(info_db[ref_gn])):
				if info_db[ref_gn][i] == []:
					info_db[ref_gn][i] = "NA"
					null_cnt += 1
				else:
					tmp_list = sorted(info_db[ref_gn][i])
					tmp_ids = []
					tmp_ids.append(tmp_list[0][-1])
					if len(tmp_list) > 1:
						for j in range(1, len(tmp_list)):
							pos, id = tmp_list[j]
							if id in tandem_db:
								tmp_ids.append(id+"-T")
							else:
								tmp_ids.append(id+"-P")
					info_db[ref_gn][i] = ','.join(tmp_ids)

			if null_cnt == allele_num:
				continue
			if is_mono.upper() == "F":
				ori_gn = ref_gn.split(":::")[0]
			else:
				ori_gn = ref_gn
			tmp_list = copy.deepcopy(ref_db[ori_gn])
			tmp_list.append(ori_gn)
			tmp_list.extend(info_db[ref_gn])
			full_allele.append(tmp_list)
		
	print("Writing allele table")
	with open(out_allele, 'w') as fout:
		allele_header = []
		for i in range(0, allele_num):
			allele_header.append("Allele %s"%(chr(i+65)))
		fout.write("#CHR\tPOS\tref gene\t%s\n"%('\t'.join(allele_header)))
		cnt = 0
		for info in sorted(full_allele):
			fout.write("%s\n"%('\t'.join(map(str, info))))
	
	print("Finished")


if __name__ == "__main__":
	if len(sys.argv) < 11:
		print("Usage: python %s <in_allele> <ref_gff3> <hap_gff3> <hap_ref_blast> <hap_hap_blast> <iden_threshold> <tandem_list> <allele_number> <out_allele> <is_mono>"%sys.argv[0])
	else:
		in_allele, ref_gff3, hap_gff3, ref_blast, hap_blast, iden, tandem, allele_num, out_allele, is_mono = sys.argv[1:]
		adjust_allele_table(in_allele, ref_gff3, hap_gff3, ref_blast, hap_blast, float(iden), tandem, int(allele_num), out_allele, is_mono)
