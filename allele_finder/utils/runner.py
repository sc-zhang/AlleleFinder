from subprocess import Popen, PIPE
from allele_finder.utils.message import Message
from allele_finder.utils.utils import FastaUtils
from allele_finder.utils.dependencies_checker import DependChecker
from os import path, makedirs, chdir, getcwd
import re


class Runner:

    def __init__(self):
        self._cmd = ""
        self._res = ""
        self._err = ""

    def set_command(self, cmd):
        self._cmd = '"%s"' % cmd

    def print_command(self):
        Message.info(self._cmd)

    def run(self):
        p = Popen(self._cmd, stdout=PIPE, stderr=PIPE, shell=True, encoding="utf-8")
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

        dep_checker = DependChecker()
        if not dep_checker.check("gmap --version"):
            Message.error("\tGMAP not found, Aborting...")
            dep_checker.exit_program_with_errors()

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
                self._cmd = "gmapl -D . -d CpDB -f 2 -n %d -t %d %s > gmap.gff3" % (
                    allele_cnt,
                    threads,
                    qry_cds,
                )
            else:
                self._cmd = "gmap -D . -d CpDB -f 2 -n %d -t %d %s > gmap.gff3" % (
                    allele_cnt,
                    threads,
                    qry_cds,
                )

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
    def blast(
            self, ref_fasta, qry_fasta, db_name, evalue, out_blast, num_alignments, threads
    ):
        if not path.exists(out_blast):
            dep_checker = DependChecker()
            if not dep_checker.check("makeblastdb -version"):
                Message.error("\tBlast not found, aborting...")
                dep_checker.exit_program_with_errors()

            Message.info("\tRunning blast")
            fa_type = FastaUtils.check_fasta_type(ref_fasta)
            self._cmd = "makeblastdb -in %s -dbtype %s -out %s" % (
                ref_fasta,
                fa_type,
                db_name,
            )
            Message.info("\tRunning command: %s" % self._cmd)
            self.run()
            with open("makeblastdb.log", "a") as fout:
                fout.write("%s\n" % self._res)
                fout.write("%s\n" % self._err)

            blast_program = "blastn" if fa_type == "nucl" else "blastp"
            self._cmd = (
                            "%s -query %s -db %s -out %s -evalue %s -outfmt 6 -num_alignments %d "
                            "-num_threads %d"
                        ) % (
                            blast_program,
                            qry_fasta,
                            db_name,
                            out_blast,
                            evalue,
                            num_alignments,
                            threads,
                        )
            Message.info("\tRunning command: %s" % self._cmd)
            self.run()
            with open("blast.log", "a") as fout:
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

        dep_checker = DependChecker()
        if not dep_checker.check("MCScanX -h"):
            Message.error("\tMCScanX not found, aborting...")
            dep_checker.exit_program_with_errors()

        if not path.exists("xyz"):
            makedirs("xyz")

        gff_id_set = set()
        if not path.exists("xyz/xyz.gff"):
            Message.info("\tLoading gff3")
            gff3_db = {}
            with open(in_gff3, "r") as fin:
                for line in fin:
                    if line.strip() == "" or line[0] == "#":
                        continue
                    data = line.strip().split()
                    if "tig" in data[0] or "ctg" in data[0]:
                        continue
                    chrn = re.findall(r"([A-Za-z]+\d+)", data[0])
                    if chrn:
                        chrn = chrn[0]
                    else:
                        continue
                    if data[2] == "gene":
                        if "Name" in data[8]:
                            gid = re.findall(r"Name=(.*)", data[8])[0].split(";")[0]
                        else:
                            gid = re.findall(r"ID=(.*)", data[8])[0].split(";")[0]

                        if chrn not in gff3_db:
                            gff3_db[chrn] = []
                        gff3_db[chrn].append([gid, int(data[3]), int(data[4])])
                        gff_id_set.add(gid)

            Message.info("\tWriting gff")
            idx = 1
            with open(path.join("xyz/xyz.gff"), "w") as fout:
                for chrn in sorted(gff3_db):
                    spid = "NN%02d" % idx
                    idx += 1
                    for gid, sp, ep in sorted(gff3_db[chrn]):
                        fout.write("%s\t%s\t%d\t%d\n" % (spid, gid, sp, ep))
        else:
            Message.info("\tGFF file found, skipping...")
            with open(path.join("xyz/xyz.gff"), "r") as fin:
                for line in fin:
                    gff_id_set.add(line.strip().split()[1])

        cds_id_set = FastaUtils.get_seq_ids(in_fa)
        match_ratio = (
                len(cds_id_set.intersection(gff_id_set))
                * 1.0
                / min(len(cds_id_set), len(gff_id_set))
        )

        if match_ratio < 0.5:
            Message.error(
                "\tFatal error: Critical low match ratio %.2f%% between gff3 file and CDS file, "
                "please check input files!" % (match_ratio * 100)
            )
            Message.error(
                "\tID in gff file: %s, ID in CDS file: %s"
                % (sorted(gff_id_set)[0], sorted(cds_id_set)[0])
            )
            Message.error(
                '\tGene ID in gff3 file are extracted with "Name" feature if exists, otherwise "ID" '
                "feature would be extracted, please check these features with CDS file"
            )
            exit(-1)

        if not path.exists("xyz/xyz.blast"):
            Message.info("\tRunning blast")
            blaster = BlastRunner()
            blaster.blast(in_fa, in_fa, "blastdb", "1e-10", "xyz/xyz.blast", 5, threads)
        else:
            Message.info("Blast file found, skipping...")

        self._cmd = "MCScanX xyz/xyz"
        Message.info("\tRunning command: %s" % self._cmd)
        self.run()

        with open("mcscanx.log", "w") as fout:
            fout.write("%s\n" % self._res)
            fout.write("%s\n" % self._err)
        Message.info("\tMCScanX finished")
        chdir(cur_dir)
