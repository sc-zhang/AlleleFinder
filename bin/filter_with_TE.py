#!/usr/bin/env python
import sys
import re


def merge_regions(regions):
	tmp_regions = []
	last_ep = 0 
	for sp, ep in sorted(regions):
		if tmp_regions == []: 
			tmp_regions.append(sp)
			last_ep = ep
		else:
			if sp > last_ep:
				tmp_regions.append(last_ep)
				tmp_regions.append(sp)
				last_ep = ep
			else:
				if ep > last_ep:
					last_ep = ep
	tmp_regions.append(last_ep)
	new_regions = []
	for i in range(0, len(tmp_regions)-1, 2): 
		new_regions.append([tmp_regions[i], tmp_regions[i+1]])
	return new_regions


def check_ovlp(id, sp, ep, pos_list, thres):	
	s = 0
	e = len(pos_list)-1
	l = -1
	while s<=e:
		mid = int((s+e)/2)
		if pos_list[mid][0] < sp:
			s = mid+1
		elif pos_list[mid][0] > sp:
			e = mid-1
		else:
			l = mid
			break
	if l==-1:
		if e==-1:
			e = 0
		if pos_list[e][1] >= sp:
			l = e
		else:
			l = e+1
	s = 0
	e = len(pos_list)-1
	r = -1
	while s<=e:
		mid = int((s+e)/2)
		if pos_list[mid][0] < ep:
			s = mid+1
		elif pos_list[mid][0] > ep:
			e = mid-1
		else:
			r = mid
			break
	if r==-1:
		if e==-1:
			e = 0
		r = e
	
	ovlp_len = 0
	if l<=r:
		tmp_regions = []
		for rsp, rep in pos_list[l: r+1]:
			ovlp = min(ep, rep)-max(sp, rsp)+1
			if ovlp <= 0:
				continue
			tmp_regions.append([max(sp, rsp), min(ep, rep)])

		if len(tmp_regions) != 0:
			tmp_regions = merge_regions(tmp_regions)
			for msp, mep in tmp_regions:
				ovlp_len += mep-msp+1
			
	if ovlp_len*1.0/(ep-sp+1) > thres:
		return True
	else:
		return False


def filter_with_TE(in_allele, hap_gff3, in_TE, TE_thres, out_allele):
	print("Loading TEs")
	TE_db = {}
	with open(in_TE, 'r') as fin:
		for line in fin:
			if line.strip() == '' or line[0] == '#':
				continue
			data = line.strip().split()
			chrn = data[0]
			sp = int(data[3])
			ep = int(data[4])
			
			if chrn not in TE_db:
				TE_db[chrn] = []
			TE_db[chrn].append([sp, ep])
	for chrn in TE_db:
		TE_db[chrn] = merge_regions(sorted(TE_db[chrn]))

	print("Getting overlap genes with TE")
	ovlp_ids = {}
	with open(hap_gff3, 'r') as fin:
		for line in fin:
			if line.strip() == '' or line[0] == '#':
				continue
			data = line.strip().split()
			if data[2] != 'gene':
				continue
			chrn = data[0]
			sp = int(data[3])
			ep = int(data[4])
			if 'Name' in data[8]:
				id = re.findall(r'Name=(.*)', data[8])[0].split(';')[0]
			else:
				id = re.findall(r'ID=(.*)', data[8])[0].split(';')[0]

			if chrn in TE_db:
				if check_ovlp(id, sp, ep, TE_db[chrn], TE_thres):
					ovlp_ids[id] = 1
	
	print("Filtering allele table")
	remove_ids = []
	with open(in_allele, 'r') as fin:
		with open(out_allele, 'w') as fout:
			for line in fin:
				if line[0] == '#':
					fout.write(line)
				else:
					data = line.strip().split()
					for i in range(3, len(data)):
						if data[i] == 'NA':
							continue
						ids = data[i].split(',')
						tmp_ids = []
						for id in ids:
							if '-' in id and id[-1] == 'P':
								rid = id.split('-')[0]
								if rid in ovlp_ids:
									remove_ids.append(rid)
									continue
							tmp_ids.append(id)
						data[i] = ','.join(tmp_ids)
					fout.write("%s\n"%('\t'.join(data)))
	remove_fn = out_allele.split('.')
	remove_fn[-1] = 'removed.txt'
	remove_fn = '.'.join(remove_fn)
	with open(remove_fn, 'w') as fout:
		fout.write('\n'.join(sorted(remove_ids)))
	print("Finished")


if __name__ == "__main__":
	if len(sys.argv) < 6:
		print("Usage: python %s <in_allele> <hap_gff3> <in_TE> <TE_threshold> <out_allele>"%sys.argv[0])
	else:
		in_allele, hap_gff3, in_TE, TE_thres, out_allele = sys.argv[1:]
		filter_with_TE(in_allele, hap_gff3, in_TE, float(TE_thres), out_allele)
		