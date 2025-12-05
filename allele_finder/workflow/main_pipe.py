import argparse
import allele_finder


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="%(prog)s {version}".format(
            version=allele_finder.__version__.__version__
        ),
    )
    subparsers = parser.add_subparsers(title="Sub commands")

    parser_construct = subparsers.add_parser("construct", help="Construct allele table")
    parser_construct.add_argument("-r", "--ref", help="reference fasta", required=True)
    parser_construct.add_argument(
        "-d", "--ref_cds", help="CDS fasta of ref", required=True
    )
    parser_construct.add_argument(
        "-f", "--ref_gff3", help="GFF3 file of ref", required=True
    )
    parser_construct.add_argument(
        "-c", "--cds", help="CDS fasta of polyploid", required=True
    )
    parser_construct.add_argument(
        "-g", "--gff3", help="GFF3 file of polyploid", required=True
    )
    parser_construct.add_argument(
        "-n", "--num_allele", help="number of allele", type=int, required=True
    )
    parser_construct.add_argument(
        "-m",
        "--is_mono",
        help="if your reference fasta is mono assembly of polyploid, add this argument",
        action="store_true",
    )
    parser_construct.add_argument(
        "--ovlp_ratio",
        help="threshold of gene pair coordinate overlap "
             "identified by GMAP, default: 0.8",
        type=float,
        default=0.8,
    )
    parser_construct.add_argument(
        "-b", "--blast_count", help="blast count, default: 2", type=int, default=2
    )
    parser_construct.add_argument(
        "--blast_identity",
        help="threshold of blast identity, default: 80",
        type=float,
        default=80,
    )
    parser_construct.add_argument(
        "-e", "--TE", help='TE gff3 for filtering, default: ""', default=""
    )
    parser_construct.add_argument(
        "-j",
        "--TE_overlap",
        help="threshold of TE overlap, default: 0.3, only effect when TE is not NULL",
        type=float,
        default=0.3,
    )
    parser_construct.add_argument(
        "--paralog_only", help="do TE filter only on paralog genes", action="store_true"
    )
    parser_construct.add_argument(
        "-w", "--workdir", help="workdir, default: wrkdir", default="wrkdir"
    )
    parser_construct.add_argument(
        "-t", "--threads", help="threads, default: 12", default=12, type=int
    )
    parser_construct.set_defaults(func=allele_finder.workflow.construct_pipe.main)

    parser_cleanup = subparsers.add_parser(
        "cleanup", help="Remove same CDS from cds, pep and gff3 files"
    )
    parser_cleanup.add_argument("--in_cds", help="Input CDS file", required=True)
    parser_cleanup.add_argument("--in_pep", help="Input PEP file")
    parser_cleanup.add_argument("--in_gff3", help="Input GFF3 file", required=True)
    parser_cleanup.add_argument("--out_cds", help="Output CDS file", required=True)
    parser_cleanup.add_argument("--out_pep", help="Output PEP file")
    parser_cleanup.add_argument("--out_gff3", help="Output GFF3 file", required=True)
    parser_cleanup.set_defaults(func=allele_finder.workflow.cleanup_pipe.main)

    parser_statistic = subparsers.add_parser("stat", help="Statistic allele table")
    parser_statistic.add_argument(
        "-i", "--input", help="Input allele table", required=True
    )
    parser_statistic.add_argument(
        "-g", "--gff3", help="GFF3 file of polyploid", required=True
    )
    parser_statistic.add_argument(
        "-o", "--output", help="Prefix of output file", required=True
    )
    parser_statistic.set_defaults(func=allele_finder.workflow.statistic_pipe.main)

    parser_adjust = subparsers.add_parser(
        "adjust", help="Adjust allele table with " "too many genes be marked as paralog"
    )
    parser_adjust.add_argument(
        "-i", "--input", help="Input allele table", required=True
    )
    parser_adjust.add_argument(
        "-m",
        "--min_num",
        help="Minium number of genes, which means the number of genes marked as paralog "
             "that distribute in different allele should be pulled down as new allele genes",
        type=int,
        required=True,
    )
    parser_adjust.add_argument(
        "-o", "--output", help="Output allele table", required=True
    )
    parser_adjust.set_defaults(func=allele_finder.workflow.adjust_pipe.main)

    try:
        args = parser.parse_args()
        args.func(args)
    except AttributeError:
        parser.print_help()
