from allele_finder.utils.runner import GmapRunner, MCScanXRunner, BlastRunner
from allele_finder.utils.message import Message
from allele_finder.utils.utils import FastaUtils, BlastUtils, AlleleUtils, GmapUtils, TEUtils
from os import path, makedirs, chdir, listdir, getcwd


def main(args):
    ref = args.ref
    ref_cds = args.ref_cds
    ref_gff3 = args.ref_gff3
    cds = args.cds
    gff3 = args.gff3
    if args.is_mono:
        is_mono = True
    else:
        is_mono = False
    allele_count = args.num_allele
    blast_count = args.blast_count
    iden_thres = args.blast_identity
    TE_file = args.TE
    TE_thres = args.TE_overlap
    wrkdir = args.workdir
    threads = args.threads
    pipeline(ref, ref_cds, ref_gff3, cds, gff3, allele_count, is_mono, blast_count, iden_thres,
             TE_file, TE_thres, wrkdir, threads)


def pipeline(ref, ref_cds, ref_gff3, cds, gff3, allele_count, is_mono, blast_count, iden_thres,
             TE_file, TE_thres, wrkdir, threads):
    if not path.exists(wrkdir):
        makedirs(wrkdir)
    ref = path.abspath(ref)
    ref_cds = path.abspath(ref_cds)
    ref_gff3 = path.abspath(ref_gff3)
    cds = path.abspath(cds)
    gff3 = path.abspath(gff3)
    if TE_file != "":
        TE_file = path.abspath(TE_file)

    Message.info("Entering: %s" % wrkdir)
    chdir(wrkdir)

    Message.info("Step1: running MCScanX")
    mcs_dir = "01.mcscanx/xyz/xyz.html"
    tandem_file = path.abspath("01.mcscanx/xyz/xyz.tandem")
    hap_blast_file = path.abspath("01.mcscanx/xyz/xyz.blast")
    if not path.exists(mcs_dir):
        mcscanxer = MCScanXRunner()
        mcscanxer.mcscanx(cds, gff3, "01.mcscanx", threads)
    else:
        Message.info("\tMCScanX result found, skipping...")

    Message.info("Step2: running GMAP")
    gmap_res = "02.gmap/gmap.gff3"
    if not path.exists(gmap_res):
        gmapper = GmapRunner()
        gmapper.gmap(ref, cds, allele_count, "02.gmap", threads)
    else:
        Message.info("\tGMAP result found, skipping...")

    Message.info("Step3: Generating first allele table")
    if not path.exists("backbone.csv"):
        Message.info("\tLoading MCScanX results")
        base_allele = []
        for fn in listdir(mcs_dir):
            full_fn = path.join(mcs_dir, fn)
            tmp_allele = AlleleUtils.get_allele_with_mcscanx(full_fn)
            base_allele.extend(tmp_allele)

        Message.info("\tLoading GMAP result")
        gff3_db, gene_order = GmapUtils.read_gff3(gmap_res)
        gff3_allele = GmapUtils.allele_gmap(gff3_db, threads)
        base_allele.extend(gff3_allele)

        Message.info("\tWriting allele list backbone")
        gff3_db, gene_order = GmapUtils.read_gff3(gff3)
        base_allele = AlleleUtils.split_allele(base_allele, gene_order)
        base_allele = AlleleUtils.merge_allele(base_allele)
        with open("backbone.csv", 'w') as fout:
            for allele in sorted(base_allele):
                fout.write("%s\n" % (",".join(sorted(allele))))
    else:
        Message.info("\tallele list backbone found, loading")
        base_allele = []
        with open("backbone.csv", 'r') as fin:
            for line in fin:
                base_allele.append(line.strip().split(','))

    backbone = path.abspath("backbone.csv")
    final_allele = base_allele

    Message.info("Step4: running blast")
    if not path.exists("03.blast"):
        makedirs("03.blast")

    curdir = getcwd()
    Message.info("\tEntering: blast")
    chdir("03.blast")

    for i in range(0, blast_count):
        Message.info("\tStarting iteration %02d" % (i + 1))
        outpre = "iter%02d" % (i + 1)
        single_fa = outpre + "_single.fa"
        multi_fa = outpre + "_multi.fa"
        out_blast = outpre + ".blast"
        if not (path.exists(single_fa) and path.exists(multi_fa)):
            FastaUtils.split_fasta_with_allele(cds, backbone, outpre)
        else:
            Message.info("\tIter %02d, Fasta file found, skipping..." % (i + 1))
        if not path.exists(out_blast):
            blaster = BlastRunner()
            blaster.blast(multi_fa, single_fa, outpre + "db", "1e-3", out_blast, 1, threads)
        else:
            Message.info("\tIter %02d, blast file found, skipping..." % (i + 1))
        final_allele.extend(BlastUtils.allele_blast(out_blast, iden_thres))
        final_allele = AlleleUtils.merge_allele(final_allele)
        backbone = outpre + ".csv"
        if not path.exists(backbone):
            with open(backbone, 'w') as fout:
                for allele in sorted(final_allele):
                    fout.write("%s\n" % (",".join(sorted(allele))))

    Message.info("\tLeaving: 03.blast")
    chdir(curdir)

    Message.info("Step5: Writing allele table")
    with open("allele.csv", 'w') as fout:
        for allele in sorted(final_allele):
            fout.write("%s\n" % (",".join(sorted(allele))))

    Message.info("Step6: Adjusting with ref annotation")
    if not (path.exists("04.ref_adjust")):
        makedirs("04.ref_adjust")

    curdir = getcwd()
    Message.info("\tEntering: 04.ref_adjust")
    chdir("04.ref_adjust")

    ref_fn = ref_cds.split('/')[-1].split('.')[0]
    hap_fn = cds.split('/')[-1].split('.')[0]
    out_blast = hap_fn + '.vs.' + ref_fn + '.blast'
    if not path.exists(out_blast):
        blaster = BlastRunner()
        blaster.blast(ref_cds, cds, "blastdb", "1e-3", out_blast, 1, threads)
    else:
        Message.info("\tBlast file found, skip")
    out_blast = path.abspath(out_blast)

    Message.info("Leaving: 04.ref_adjust")
    chdir(curdir)

    allele_file = "allele.adjusted.txt"
    if not path.exists(allele_file):
        AlleleUtils.adjust_allele_table("allele.csv", ref_gff3, gff3, out_blast, hap_blast_file, iden_thres,
                                        tandem_file, allele_count, allele_file, is_mono, "04.ref_adjust.log")
    else:
        Message.info("\tallele.adjusted.txt found, skipping...")

    if TE_file != "":
        Message.info("Step7: Filtering with TEs")
        allele_file = "allele.adjusted.nonTEs.txt"
        if not path.exists(allele_file):
            TEUtils.filter_with_TE("allele.adjusted.txt", gff3, TE_file, TE_thres,
                                   "allele.adjusted.nonTEs.txt", "05.filter.log")
        else:
            Message.info("\tallele.adjusted.nonTEs.txt found, skipping...")

    Message.info("Allele table constructed")
