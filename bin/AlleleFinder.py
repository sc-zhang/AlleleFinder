#!/usr/bin/env python
import sys
import os
import argparse
import time
import allele_backbone as ab
import allele_gmap as ag
import allele_blast as abl


def time_print(info):
	print("\033[32m%s\033[0m %s"%(time.strftime('[%H:%M:%S]',time.localtime(time.time())), info))


def get_opts():
	group = argparse.ArgumentParser()
	group.add_argument('-m', '--mono', help="reference mono fasta", required=True)
	group.add_argument('-c', '--cds', help="CDS fasta", required=True)
	group.add_argument('-n', '--num_allele', help="number of allele", type=int, required=True)
	group.add_argument('-g', '--gff3', help="GFF3 file", required=True)
	group.add_argument('-w', '--workdir', help="workdir, default: wrkdir", default="wrkdir")
	group.add_argument('-t', '--threads', help="threads, default: 12", default=12, type=int)
	return group.parse_args()


def AlleleFinder(mono, cds, gff3, na, wrkdir, threads):
	if not os.path.exists(wrkdir):
		os.mkdir(wrkdir)
	
	mono = os.path.abspath(mono)
	cds = os.path.abspath(cds)
	gff3 = os.path.abspath(gff3)

	script_dir = sys.path[0]
	time_print("Entering: %s"%wrkdir)
	os.chdir(wrkdir)

	time_print("Step1: running MCScanX")
	mcs_dir = "mcscanx/xyz/xyz.html"
	if not os.path.exists(mcs_dir):
		cmd = "%s/run_MCScanX.py %s %s mcscanx %d &> step1.log"%(script_dir, cds, gff3, threads)
		time_print("\tRunning: %s"%cmd)
		os.system(cmd)
	else:
		time_print("\tMCScanX result found, skip")
	
	time_print("Step2: running gmap")
	gmap_res = "gmap/gmap.gff3"
	if not os.path.exists(gmap_res):
		cmd = "%s/run_gmap.py %s %s %d gmap %d &> step2.log"%(script_dir, mono, cds, na, threads)
		time_print("\tRunning: %s"%cmd)
		os.system(cmd)
	else:
		time_print("\tGmap result found, skip")
	
	time_print("Step3: Generating allele table")
	time_print("\tLoading MCScanX results")
	base_allele = []
	for fn in os.listdir(mcs_dir):
		full_fn = os.path.join(mcs_dir, fn)
		time_print("\tLoading: %s"%full_fn)
		tmp_allele = ab.get_allele_with_mcscanx(full_fn)
		base_allele.extend(tmp_allele)
	
	time_print("\tLoading GMAP result")
	gff3_db, gene_order = ag.read_gff3(gmap_res)
	gff3_allele = ag.allele_gmap(gff3_db, threads)
	base_allele.extend(gff3_allele)

	time_print("\tWriting allele list")
	gff3_db, gene_order = ag.read_gff3(gff3)
	base_allele = ab.split_allele(base_allele, gene_order)
	allele_list = ab.merge_allele(base_allele)
	with open("result.csv", 'w') as fout:
		for allele in sorted(allele_list):
			fout.write("%s\n"%(",".join(sorted(allele))))
	
	time_print("Finished")


if __name__ == "__main__":
	opts = get_opts()
	mono = opts.mono
	cds = opts.cds
	na = opts.num_allele
	gff3 = opts.gff3
	wrkdir = opts.workdir
	threads = opts.threads
	AlleleFinder(mono, cds, gff3, na, wrkdir, threads)
