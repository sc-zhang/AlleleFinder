#!/usr/bin/env python
import sys
import os


def run_gmap(ref_mono, in_cds, allele_cnt, workdir, threads):
	ref_mono = os.path.abspath(ref_mono)
	in_cds = os.path.abspath(in_cds)

	if not os.path.exists(workdir):
		os.mkdir(workdir)
	print("Entering %s"%workdir)
	os.chdir(workdir)
	
	if not os.path.exists("CpDB"):
		cmd = "gmap_build -D . -d CpDB %s &> /dev/null"%ref_mono
		print("Running command: %s"%cmd)
		os.system(cmd)
	else:
		print("Gmap db found, skip")
	
	if not os.path.exists("gmap.gff3"):
		fs = os.path.getsize(ref_mono)
		if fs >= 2**32:
			cmd = "gmapl -D . -d CpDB -f 2 -n %s -t %s %s 1> gmap.gff3 2>/dev/null"%(allele_cnt, threads, in_cds)
		else:
			cmd = "gmap -D . -d CpDB -f 2 -n %s -t %s %s 1> gmap.gff3 2>/dev/null"%(allele_cnt, threads, in_cds)
		
		
		print("Running command: %s"%cmd)
		os.system(cmd)
	else:
		print("gmap result found, skip")
	
	print("Finished")


if __name__ == "__main__":
	if len(sys.argv) < 6:
		print("Usage: python %s <ref_mono> <in_cds> <allele_count> <workdir> <threads>"%sys.argv[0])
	else:
		ref_mono, in_cds, allele_cnt, workdir, threads = sys.argv[1:]
		run_gmap(ref_mono, in_cds, allele_cnt, workdir, threads)
