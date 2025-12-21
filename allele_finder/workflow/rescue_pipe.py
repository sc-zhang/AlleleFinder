from allele_finder.utils.utils import FastaUtils
from allele_finder.utils.message import Message
import re


def rescue(in_allele_file, in_gff3_file, in_cds_file, in_pep_file, out_allele_file):
    Message.info("Loading allele genes")
    allele_table = []
    exists_allele_genes = set()
    allele_table_header = ""
    with open(in_allele_file, "r") as fin:
        for line in fin:
            if line[0] == "#":
                allele_table_header = line.strip()
            else:
                data = line.strip().split()
                allele_table.append(data)
                for allele in data[3:]:
                    if allele == "NA":
                        continue
                    for sub_gene in allele.split(","):
                        if "-" in sub_gene:
                            gid = sub_gene.split("-")[0]
                        else:
                            gid = sub_gene
                        exists_allele_genes.add(gid)

    Message.info("Loading hap gff3")
    hap_db = {}
    with open(in_gff3_file, "r") as fin:
        for line in fin:
            if line.strip() == "" or line[0] == "#":
                continue
            data = line.strip().split()
            if not data[0].upper().startswith("CHR"):
                continue
            chrn = data[0][:-1]
            hap = ord(data[0][-1]) - 65
            try:
                if "Name" in data[8]:
                    gid = re.findall(r"Name=(.*)", data[8])[0].split(";")[0]
                else:
                    gid = re.findall(r"ID=(.*)", data[8])[0].split(";")[0]
                hap_db[gid] = chrn, hap
            except IndexError:
                continue

        if not hap_db:
            Message.error("No ID or Name found in %s" % in_gff3_file)
            exit(-1)

    Message.info("Getting genes with same CDS sequences")
    # if genes with full same CDS sequences, mark as identical genes
    identical_gene_set, identical_gene_db = FastaUtils.get_full_same_seqs(in_cds_file, True)
    identical_group_matched_exists_genes = {}
    for gene in exists_allele_genes:
        if gene in identical_gene_set:
            for group in identical_gene_db:
                if gene in identical_gene_db[group]:
                    identical_group_matched_exists_genes[gene] = group
                    break

    partial_identical_gene_set = set()
    partial_identical_gene_db = {}
    partial_identical_group_matched_exists_genes = {}

    if in_pep_file:
        # if genes with same PEP sequences but had different base with CDS sequences, mark as allele genes
        Message.info("Getting genes with same PEP sequences")
        partial_identical_gene_set, partial_identical_gene_db = FastaUtils.get_full_same_seqs(
            in_pep_file, True
        )
        for gene in exists_allele_genes:
            if gene in partial_identical_gene_set:
                for group in partial_identical_gene_db:
                    if gene in partial_identical_gene_db[group]:
                        partial_identical_group_matched_exists_genes[gene] = group
                        break

    Message.info("Adding cleaned genes and writing new allele table")
    with open(out_allele_file, "w") as fout:
        fout.write("%s\n" % allele_table_header)
        for row in range(len(allele_table)):
            sub_genes = set()
            cur_row = []
            for allele in allele_table[row][3:]:
                cur_row.append(allele.split(","))
                if allele == "NA":
                    continue
                for sub_gene in allele.split(","):
                    if "-" in sub_gene:
                        gid = sub_gene.split("-")[0]
                    else:
                        gid = sub_gene
                    sub_genes.add(sub_gene)

            # get identical genes and allele genes which could be added to current allele table
            cur_identical_genes = set()
            cur_partial_identical_genes = set()
            for gene in sub_genes:
                if gene in identical_group_matched_exists_genes:
                    group = identical_group_matched_exists_genes[gene]
                    cur_identical_genes.update(identical_gene_db[group])
                if gene in partial_identical_group_matched_exists_genes:
                    group = partial_identical_group_matched_exists_genes[gene]
                    for gid in partial_identical_gene_db[group]:
                        if gid not in identical_gene_set:
                            cur_partial_identical_genes.add(gid)

            # add identical genes
            for gene in cur_identical_genes:
                if gene not in hap_db:
                    continue
                _, hap = hap_db[gene]
                if cur_row[hap][0] == "NA":
                    cur_row[hap][0] = gene
                else:
                    cur_row[hap].append(gene + "-P")

            # add allele genes
            for gene in cur_partial_identical_genes:
                if gene not in hap_db:
                    continue
                _, hap = hap_db[gene]
                if cur_row[hap][0] == "NA":
                    cur_row[hap][0] = gene
                else:
                    cur_row[hap].append(gene + "-P")

            info = allele_table[row][:3]
            for _ in cur_row:
                info.append(",".join(_))
            fout.write("%s\n" % ("\t".join(info)))

    Message.info("Finished")


def main(args):
    in_allele = args.input
    in_gff3 = args.gff3
    in_cds = args.cds
    in_pep = args.pep
    out_allele = args.output
    rescue(in_allele, in_gff3, in_cds, in_pep, out_allele)
