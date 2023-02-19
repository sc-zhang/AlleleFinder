#!/usr/bin/env python
import sys
import re
import copy
import dup_pipe as dp

def adjust_allele_table(in_allele, ref_gff3, hap_gff3, ref_blast, hap_blast, iden, tandem, allele_num, out_allele, is_mono, pFile, tFile):
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
	line_cnt = 0
	with open(in_allele, 'r') as fin:
		for line in fin:
			data = line.strip().split(',')
			line_cnt += 1
			for id in data:
				hchrn, hpos = hap_db[id]
				hidx = ord(hchrn[-1])-65
				if id not in ref_blast_db:
					if id not in hap_blast_db:
						continue  # 这种一般是singleton基因
					ref_id = hap_blast_db[id][0]
				else:
					ref_id = id
				if ref_id not in ref_blast_db or ref_id not in data:  # id不在等位中或者比对不上reference
					continue
				ref_gn = ref_blast_db[ref_id][0]
				if ref_gn not in ref_db:
					continue
				if is_mono.upper() == 'T':
					if ref_gn not in info_db:
						info_db[ref_gn] = [[] for i in range(0, allele_num)]
					info_db[ref_gn][hidx].append([hpos, id])
				else:
					tmp_gn = "%s:::%d"%(ref_gn, line_cnt)
					if tmp_gn not in info_db:
						info_db[tmp_gn] = [[] for i in range(0, allele_num)]
					info_db[tmp_gn][hidx].append([hpos, id])
		
		dup_gene_db = {}
		allele2ref = {}
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
							dup_gene_db[id]
						tmp_ids = tmp_list[0][-1]
						
#							if id in tandem_db:
#								tmp_ids.append(id+"-T")
#							else:
#								tmp_ids.append(id+"-P")
					#allele2ref[tmp_ids] = ref_g
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
			for tmp_g in info_db[ref_gn]:
				allele2ref[tmp_g] = ref_gn
			full_allele.append(tmp_list)
	#### 开始重新分配被拉取下来的基因	
		singleton = {}
		allele2table = {}
		addition_allele = {}
		dup_classter = {}
		network2gn = {}
		gn2network = {}
		nc = 0
		for g in dup_gene_db:
			if g not in hap_blast_db:
				singleton[g] = ''
				continue
			else:
				blast_gn = hap_blast_db[g][0]
				if blast_gn in dup_gene_db:
					if blast_gn not in gn2network and g not in gn2network:
						nc += 1
						network2gn[nc] = set([blast_gn, g])
						gn2network[blast_gn] = nc
						gn2network[g] = nc
					elif g in gn2network:
						nc = gn2network[g]
						network2gn[nc].add(blast_gn)
					else:
						nc = gn2network[blast_gn]
						network2gn[nc].add(g)
				elif blast_gn in allele2ref:
					ref_gn = allele2ref[blast_gn]
					hchrn, hpos = hap_db[id]
					hidx = ord(hchrn[-1])-65
					candidate_gn = info_db[ref_gn][hidx]
					dup_classter[g] = [ref_gn, 'remove']
					if candidate_gn == 'NA':
						info_db[ref_gn][hidx] = g
						dup_classter[g][1] = 'keep'
				else:
					continue
		gnNet = list(network2gn.values())
		for net in gnNet:
			# check if g in allele table
			tmp_ref_gn = 'NA'
			for tmpg in net:
				if tmpg in dup_classter:
					tmp_ref_gn = dup_classter[g][0]
			if tmp_ref_gn == 'NA':
				ref_gn = list(net)[0]
				info_db[ref_gn] = ['NA' for i in range(0, allele_num)]
				for tmpg in net:
					hchrn, hpos = hap_db[tmpg]
					hidx = ord(hchrn[-1])-65
					info_db[ref_gn][hidx] = tmpg
				addition_allele[ref_gn] = ''
			else:
				tmpAddition = []
				for tmpg in net: # check every gene
					hchrn, hpos = hap_db[tmpg]
					hidx = ord(hchrn[-1])-65
					candidate_gn = info_db[tmp_ref_gn][hidx]
					# 其中几个基因能被拉回等位表中
					if candidate_gn == 'NA':
						info_db[tmp_ref_gn][hidx] = tmpg
						if tmp_ref_gn not in allele2table:
							allele2table[tmp_ref_gn] = ''
					elif dup_classter[tmpg][1] == 'keep':
						# 避免能挂到等位表上的基因被去掉
						if tmp_ref_gn not in allele2table:
							allele2table[tmp_ref_gn] = ''
					else:
						tmpAddition.append(tmpg)
				if len(tmpAddition) >= 1:
					ref_gn = tmpAddition[0]
					info_db[ref_gn] = ['NA' for i in range(0, allele_num)]
				for tmpg in tmpAddition:
					hchrn, hpos = hap_db[tmpg]
					hidx = ord(hchrn[-1])-65
					info_db[ref_gn][hidx] = tmpg
				addition_allele[ref_gn] = ''
				# 将挂不上等位表的基因自成一行
			
		for dup_in_allele in dup_classter:
			ref_gn = dup_classter[dup_in_allele][0]
			if ref_gn not in allele2table and dup_in_allele not in gn2network:
				hchrn, hpos = hap_db[dup_in_allele]
				hidx = ord(hchrn[-1])-65
				info_db[ref_gn] = ['NA' for i in range(0, allele_num)]
				candidate_gn = info_db[dup_in_allele][hidx]
				if candidate_gn == 'NA':
					info_db[ref_gn][hidx] = dup_in_allele
					allele2table[ref_gn] = ''
				else:
					ref_gn = dup_in_allele
					info_db[ref_gn][hidx] = dup_in_allele
					addition_allele[ref_gn] = ''
		for sg in singleton:
			ref_gn = sg
			hchrn, hpos = hap_db[sg]
			hidx = ord(hchrn[-1])-65
			info_db[ref_gn] = ['NA' for i in range(0, allele_num)]
			info_db[ref_gn][hidx] = sg
	### 重新加上paralog的基因
	### load pDic & tDic
	pDic = dp.read_repeat_file(pFile)
	tDic = dp.read_repeat_file(tFile)
	full_allele = dp.add_dup_info(full_allele, info_db, hap_db, singleton, allele2table, addition_allele, pDic, tDic)
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
		print("Usage: python %s <in_allele> <ref_gff3> <hap_gff3> <hap_ref_blast> <hap_hap_blast> <iden_threshold> <tandem_list> <allele_number> <out_allele> <is_mono> <pFile> <tFile>"%sys.argv[0])
	else:
		in_allele, ref_gff3, hap_gff3, ref_blast, hap_blast, iden, tandem, allele_num, out_allele, is_mono, pFile, tFile = sys.argv[1:]
		adjust_allele_table(in_allele, ref_gff3, hap_gff3, ref_blast, hap_blast, float(iden), tandem, int(allele_num), out_allele, is_mono, pFile, tFile)
