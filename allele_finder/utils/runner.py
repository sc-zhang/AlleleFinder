from subprocess import Popen, PIPE
from allele_finder.utils.message import Message
from allele_finder.utils.utils import FastaUtils
from os import path, makedirs, chdir, getcwd
import re


class Runner:

    def __init__(self):
        self._cmd = ""
        self._res = ""
        self._err = ""

    def set_command(self, cmd):
        self._cmd = "\"%s\"" % cmd

    def print_command(self):
        Message.info(self._cmd)

    def run(self):
        p = Popen(self._cmd, stdout=PIPE, stderr=PIPE, shell=True, encoding='utf-8')
        self._res, self._err = p.communicate()

    def get_result(self):
        return self._res

    def get_err(self):
        return self._err


class GmapRunner(Runner):
    def gmap(self, ref_mono, qry_cds, allele_cnt, workdir, threads):
        ref_mono = path.abspath(ref_mono)
        qry_cds = path.abspath(qry_cds)
        cur_dir = getcwd()
        if not path.exists(workdir):
            makedirs(workdir)
        Message.info("\tEntering %s" % workdir)
        chdir(workdir)

        if not path.exists("CpDB"):
            self._cmd = "gmap_build -D . -d CpDB %s" % ref_mono
            Message.info("\tRunning command: %s" % self._cmd)
            self.run()
            with open("gmap_build.log", "w") as fout:
                fout.write("%s\n" % self._res)
                fout.write("%s\n" % self._err)
        else:
            Message.info("\tGmap db found, skipping...")

        if not path.exists("gmap.gff3"):
            fs = path.getsize(ref_mono)
            if fs >= 2 ** 32:
                self._cmd = "gmapl -D . -d CpDB -f 2 -n %d -t %d %s > gmap.gff3" % (allele_cnt, threads, qry_cds)
            else:
                self._cmd = "gmap -D . -d CpDB -f 2 -n %d -t %d %s > gmap.gff3" % (allele_cnt, threads, qry_cds)

            Message.info("\tRunning command: %s" % self._cmd)
            self.run()
            with open("gmap.log", "w") as fout:
                fout.write("%s\n" % self._res)
                fout.write("%s\n" % self._err)
        else:
            Message.info("\tgmap result found, skipping...")

        Message.info("\tGmap finished")
        chdir(cur_dir)


class BlastRunner(Runner):
    def blast(self, ref_fasta, qry_fasta, db_name, out_blast, threads):
        if not path.exists(out_blast):
            Message.info("\tRunning blast")
            fa_type = FastaUtils.check_fasta_type(ref_fasta)
            self._cmd = "makeblastdb -in %s -dbtype %s -out %s" % (ref_fasta, fa_type, db_name)
            Message.info("\tRunning command: %s" % self._cmd)
            self.run()
            with open("makeblastdb.log", 'a') as fout:
                fout.write("%s\n" % self._res)
                fout.write("%s\n" % self._err)

            self._cmd = ("blastn -query %s -db %s -out %s -evalue 1e-10 -outfmt 6 -num_alignments 5 "
                         "-num_threads %d") % (qry_fasta, db_name, out_blast, threads)
            Message.info("\tRunning command: %s" % self._cmd)
            self.run()
            with open("blast.log", 'a') as fout:
                fout.write("%s\n" % self._res)
                fout.write("%s\n" % self._err)
        else:
            Message.info("\tBlast file found, skipping...")

        Message.info("\tBlast finished")


class MCScanXRunner(Runner):
    def mcscanx(self, in_fa, in_gff3, out_dir, threads):
        in_fa = path.abspath(in_fa)
        in_gff3 = path.abspath(in_gff3)
        if not path.exists(out_dir):
            makedirs(out_dir)

        cur_dir = getcwd()
        Message.info("\tEntering %s" % out_dir)
        chdir(out_dir)
        if not path.exists("xyz"):
            makedirs("xyz")

        if not path.exists("xyz/xyz.blast"):
            Message.info("\tRunning blast")
            blaster = BlastRunner()
            blaster.blast(in_fa, in_fa, "blastdb", "xyz/xyz.blast", threads)
        else:
            Message.info("Blast file found, skipping...")

        if not path.exists("xyz/xyz.gff"):
            Message.info("\tLoading gff3")
            gff3_db = {}
            with open(in_gff3, 'r') as fin:
                for line in fin:
                    if line.strip() == '' or line[0] == '#':
                        continue
                    data = line.strip().split()
                    chrn = re.findall(r'([A-Za-z]+\d+)', data[0])[0]
                    if 'tig' in chrn or 'ctg' in chrn:
                        continue
                    if data[2] == 'gene':
                        if "Name" in data[8]:
                            gid = re.findall(r'Name=(.*)', data[8])[0].split(';')[0]
                        else:
                            gid = re.findall(r'ID=(.*)', data[8])[0].split(';')[0]

                        if chrn not in gff3_db:
                            gff3_db[chrn] = []
                        gff3_db[chrn].append([gid, int(data[3]), int(data[4])])

            Message.info("\tWriting gff")
            idx = 1
            with open(path.join("xyz/xyz.gff"), 'w') as fout:
                for chrn in sorted(gff3_db):
                    spid = "NN%02d" % (idx)
                    idx += 1
                    for gid, sp, ep in sorted(gff3_db[chrn]):
                        fout.write("%s\t%s\t%d\t%d\n" % (spid, gid, sp, ep))
        else:
            Message.info("\tGFF file found, skipping...")

        self._cmd = "MCScanX xyz/xyz"
        Message.info("\tRunning command: %s" % self._cmd)
        self.run()

        with open("mcscanx.log", 'w') as fout:
            fout.write("%s\n" % self._res)
            fout.write("%s\n" % self._err)
        Message.info("\tMCScanX finished")
        chdir(cur_dir)
