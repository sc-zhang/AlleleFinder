import argparse
import allele_finder


def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(title="Sub commands")
    parser_prepare = subparsers.add_parser('prepare', help="Remove same CDS from cds, pep and gff3 files")
    parser_prepare.add_argument('--in_cds', help="Input CDS file", required=True)
    parser_prepare.add_argument('--in_pep', help="Input PEP file")
    parser_prepare.add_argument('--in_gff3', help="Input GFF3 file", required=True)
    parser_prepare.add_argument('--out_cds', help="Output CDS file", required=True)
    parser_prepare.add_argument('--out_pep', help="Output PEP file")
    parser_prepare.add_argument('--out_gff3', help="Output GFF3 file", required=True)
    parser_prepare.set_defaults(func=allele_finder.workflow.prepare_pipe.main)

    parser_construct = subparsers.add_parser('construct', help="Construct allele table")
    parser_construct.add_argument('-r', '--ref', help="reference fasta", required=True)
    parser_construct.add_argument('-d', '--ref_cds', help="CDS fasta of ref", required=True)
    parser_construct.add_argument('-f', '--ref_gff3', help="GFF3 file of ref", required=True)
    parser_construct.add_argument('-c', '--cds', help="CDS fasta of polyploid", required=True)
    parser_construct.add_argument('-g', '--gff3', help="GFF3 file of polyploid", required=True)
    parser_construct.add_argument('-n', '--num_allele', help="number of allele", type=int, required=True)
    parser_construct.add_argument('-m', '--is_mono',
                                  help="If your reference fasta is mono assembly of polyploid, add this argument",
                                  action="store_true")
    parser_construct.add_argument('-b', '--blast_count', help="blast count, default: 2",
                                  type=int, default=2)
    parser_construct.add_argument('-i', '--blast_identity', help="threshold of blast identity, "
                                                                 "default: 80", type=float, default=80)
    parser_construct.add_argument('-e', '--TE', help="TE gff3 for filtering, default: \"\"", default="")
    parser_construct.add_argument('-j', '--TE_overlap',
                                  help="threshold of TE overlap, default: 0.3, only effect when TE is not NULL",
                                  type=float,
                                  default=0.3)
    parser_construct.add_argument('-w', '--workdir', help="workdir, default: wrkdir", default="wrkdir")
    parser_construct.add_argument('-t', '--threads', help="threads, default: 12", default=12, type=int)
    parser_construct.set_defaults(func=allele_finder.workflow.construct_pipe.main)

    parser_statistic = subparsers.add_parser('stat', help="Statistic allele table")
    parser_statistic.add_argument('-i', '--input', help="Input allele table", required=True)
    parser_statistic.add_argument('-g', '--gff3', help="GFF3 file of polyploid", required=True)
    parser_statistic.add_argument('-o', '--output', help="Prefix of output file", required=True)
    parser_statistic.set_defaults(func=allele_finder.workflow.statistic_pipe.main)

    parser_adjust = subparsers.add_parser('adjust', help="Adjust allele table with "
                                                         "too many genes be marked as paralog")
    parser_adjust.add_argument('-i', '--input', help="Input allele table", required=True)
    parser_adjust.add_argument('-m', '--min_num',
                               help="Minium number of genes, which means the number of genes marked as paralog "
                                    "that distribute in different allele should be pulled down as new allele genes",
                               type=int, required=True)
    parser_adjust.add_argument('-o', '--output', help="Output allele table", required=True)
    parser_adjust.set_defaults(func=allele_finder.workflow.adjust_pipe.main)

    try:
        args = parser.parse_args()
        args.func(args)
    except AttributeError as e:
        print(e)
        parser.print_help()
