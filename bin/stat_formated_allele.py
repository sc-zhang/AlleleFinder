#!/usr/bin/env python
import sys
import re


def stat_formated_allele(in_allele, in_hap_gff3, out_stat):
	print("Loading hap gff3")
	hap_db = {}
	with open(in_hap_gff3, 'r') as fin:
		for line in fin:
			if line.strip() == '' or line[0] == '#':
				continue
			data = line.strip().split()
			if 'ctg' in data[0] or 'tig' in data[0] or 'utg' in data[0]:
				continue
			chrn = "Chr%02d"%(int(re.findall(r'[A-Za-z](\d+)', data[0])[0]))
			if "Name" in data[8]:
				id = re.findall(r'Name=(.*)', data[8])[0].split(';')[0]
			else:
				id = re.findall(r'ID=(.*)', data[8])[0].split(';')[0]
			hap_db[id] = chrn
	
	print("Loading allele table")
	allele_summary = {}
	with open(in_allele, 'r') as fin:
		for line in fin:
			data = line.strip().split('\t')
			if line[0] == '#':
				allele_cnt = len(data)-3
				continue

			tand_cnt = 0
			disp_cnt = 0
			allele_cnt_list = [0 for i in range(0, allele_cnt)]
			for i in range(3, len(data)):
				if data[i] == 'NA':
					continue
				ids = data[i].split(',')
				chrn_cnt_db = {}
				allele_cnt_list[i-3] += 1
				for id in ids:
					if '-' in id:
						id, type = id.split('-')
						if type == 'T':
							tand_cnt += 1
						else:
							disp_cnt += 1
					if id in hap_db:
						chrn = hap_db[id]
						if chrn not in chrn_cnt_db:
							chrn_cnt_db[chrn] = 0
						chrn_cnt_db[chrn] += 1
			max_chrn = ""
			max_cnt = 0
			for chrn in chrn_cnt_db:
				if chrn_cnt_db[chrn] > max_cnt:
					max_cnt = chrn_cnt_db[chrn]
					max_chrn = chrn
			if max_chrn not in allele_summary:
				# [[allele count list], tandem count, dispersely count, line count]
				allele_summary[max_chrn] = [[0 for i in range(0, allele_cnt)], 0, 0, 0]
			
			allele_summary[max_chrn][-1] += 1
			for i in range(0, allele_cnt):
				allele_summary[max_chrn][0][i] += allele_cnt_list[i]
			
			allele_summary[max_chrn][1] += tand_cnt
			allele_summary[max_chrn][2] += disp_cnt
	
	print("Writing summary")
	with open(out_stat, 'w') as fout:
		allele_list = []
		for i in range(0, allele_cnt):
			allele_list.append("Allele %s"%chr(65+i))
		fout.write("#CHR,Line count,%s,Tandem,Dispersely\n"%(','.join(allele_list)))
		sum_info = [0 for i in range(0, allele_cnt+3)]
		for chrn in sorted(allele_summary):
			info = [chrn, allele_summary[chrn][3]]
			info.extend(allele_summary[chrn][0])
			info.append(allele_summary[chrn][1])
			info.append(allele_summary[chrn][2])
			sum_info[0] += allele_summary[chrn][3]
			for i in range(0, len(allele_summary[chrn][0])):
				sum_info[i+1] += allele_summary[chrn][0][i]
			sum_info[allele_cnt+1] += allele_summary[chrn][1]
			sum_info[allele_cnt+2] += allele_summary[chrn][2]
			fout.write("%s\n"%(','.join(map(str, info))))
		fout.write("Total,%s\n"%(','.join(map(str, sum_info))))
	print("Finished")


if __name__ == "__main__":
	if len(sys.argv) < 4:
		print("Usage: python %s <in_allele> <in_hap_gff3> <out_stat>"%sys.argv[0])
	else:
		in_allele, in_hap_gff3, out_stat = sys.argv[1:]
		stat_formated_allele(in_allele, in_hap_gff3, out_stat)
