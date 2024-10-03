from copy import deepcopy
from allele_finder.utils.message import Message


def pull_down_paralog(in_allele, min_num, out_allele):
    Message.info("Pulling down")
    header_lines = []
    allele_data = []
    with open(in_allele, 'r') as fin:
        for line in fin:
            # Get allele counts, the header of allele table is like below:
            # #CHR  POS ref gene    Allele A    Allele B
            # different columns are divided by tab
            if line[0] == '#':
                data = line.strip().split('\t')
                allele_cnt = len(data[3:])
                header_lines.append(line)
            else:
                data = line.strip().split('\t')
                allele_data.append(data)

    is_complete = False
    round_idx = 1
    while not is_complete:
        Message.info("Round %d" % round_idx)
        round_idx += 1
        new_allele_data = []
        is_complete = True
        for data in allele_data:
            # Count paralogs distribute in different alleles
            para_in_diff_allele = 0
            for i in range(3, len(data)):
                genes = data[i]
                # Each column may contain more than one gene, and the first one is allelic gene,
                # and the others may be marked with "-T" as tandem and "-P" as paralog
                sub_genes = genes.split(',')
                for gene in sub_genes:
                    # Check each allele has paralog gene or not, every allele only count once
                    if gene.split('-')[-1] == 'P':
                        para_in_diff_allele += 1
                        break

            # if the number of alleles which have paralog genes greater than threshold,
            # the first paralog gene in each allele should be pulled down as new alleles
            new_allele = []
            if para_in_diff_allele >= min_num:
                is_complete = False
                new_allele = data[:3]
                for _ in range(allele_cnt):
                    new_allele.append("NA")
                for i in range(3, len(data)):
                    sub_genes = data[i].split(',')
                    new_sub_genes = []
                    # Only pull down first paralog
                    is_fst = True
                    for j in range(len(sub_genes)):
                        if is_fst and sub_genes[j].split('-')[-1] == 'P':
                            new_allele[i] = '-'.join(sub_genes[j].split('-')[:-1])
                            is_fst = False
                        else:
                            new_sub_genes.append(sub_genes[j])
                    # Update current allele
                    data[i] = ','.join(new_sub_genes)

            if new_allele:
                new_allele_data.append(new_allele)
            new_allele_data.append(data)
        allele_data = deepcopy(new_allele_data)

    with open(out_allele, 'w') as fout:
        fout.write("%s" % ''.join(header_lines))
        for data in allele_data:
            fout.write("%s\n" % ('\t'.join(data)))

    Message.info("Finished")


def main(args):
    in_allele = args.input
    min_num = args.min_num
    out_allele = args.output
    pull_down_paralog(in_allele, min_num, out_allele)
