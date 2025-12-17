from allele_finder.utils.utils import FastaUtils, GFF3Utils
from allele_finder.utils.message import Message


def cleanup_files(in_cds_file, in_pep_file, in_gff3_file, out_cds_file, out_pep_file, out_gff3_file, by_pep):
    if not in_cds_file and not in_pep_file:
        Message.error("At least one of CDS or PEP file must be provided.")
        exit(-1)
    if not in_pep_file and by_pep:
        Message.error("Could not cleanup by pep while no PEP file provided.")
        exit(-1)

    Message.info("Getting duplicate Sequences")
    if in_cds_file and not by_pep:
        retain_db, dup_db = FastaUtils.get_full_same_seqs(in_cds_file)
    else:
        retain_db, dup_db = FastaUtils.get_full_same_seqs(in_pep_file)

    drop_set = set()
    for gid in retain_db:
        if retain_db[gid] != gid:
            drop_set.add(gid)
    Message.info("%d duplicate sequences found" % len(drop_set))

    with open("DuplicateSeq.list", 'w') as fout:
        fout.write("#RetainID\tDupID\n")
        for gid in sorted(dup_db):
            fout.write("%s\t%s\n" % (gid, ','.join(dup_db[gid])))

    Message.info("Filtering CDS file")
    FastaUtils.filter_fasta(in_cds_file, out_cds_file, drop_set)

    if in_pep_file:
        Message.info("Filtering PEP file")
        FastaUtils.filter_fasta(in_pep_file, out_pep_file, drop_set)

    Message.info("Filtering GFF3 file")
    GFF3Utils.filter_gff3(in_gff3_file, out_gff3_file, drop_set)

    Message.info("Finished")


def main(args):
    in_cds_file = args.in_cds
    in_pep_file = args.in_pep
    in_gff3_file = args.in_gff3
    out_cds_file = args.out_cds
    out_pep_file = args.out_pep
    out_gff3_file = args.out_gff3
    by_pep = args.by_pep
    cleanup_files(in_cds_file, in_pep_file, in_gff3_file, out_cds_file, out_pep_file, out_gff3_file, by_pep)
