#!/usr/bin/env python
import argparse


def get_opts():
    group = argparse.ArgumentParser()
    group.add_argument('-i', '--input', help="Input allele table", required=True)
    group.add_argument('-m', '--min_num', help="Minium number of genes, which means the number of genes marked as paralog that distribute in different allele should be pulled down as new allele genes", type=int, required=True)
    group.add_argument('-o', '--output', help="Output allele table", required=True)


def pull_down_paralog(in_allele, min_num, out_allele):
    print("Pulling down")
    with open(in_allele, 'r') as fin:
        with open(out_allele, 'w') as fout:
            for line in fin:
                # Get allele counts, the header of allele table is like below:
                # #CHR  POS ref gene    Allele A    Allele B
                # different colunms are divided by tab
                if line[0] == '#':
                    data = line.strip().split('\t')
                    allele_cnt = len(data[3:])
                    # Copy the header line
                    fout.write(line)
                
                # Start dealing with allele table
                else:
                    data = line.strip().split('\t')
                    # Count paralogs distribute in different alleles
                    para_in_diff_allele = 0
                    for i in range(3, len(data)):
                        genes = data[i]
                        # Each column may contain more than one gene, and the first one is allelic gene, and the others may marked with "-T" as tandem and "-P" as paralog
                        sub_genes = genes.split(',')
                        for gene in sub_genes:
                            # Check each allele has paralog gene or not, every allele only count once
                            if gene.split('-')[-1] == 'P':
                                para_in_diff_allele += 1
                                break
                        
                    # if the number of alleles which have paralog genes greater than threshold, the first paralog gene in each allele should be pulled down as new alleles
                    if para_in_diff_allele >= min_num:
                        new_allele = data[:3]
                        for _ in range(allele_cnt):
                            new_allele.append("")
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
                        
                        # Write alleles which after pulled down
                        fout.write("%s\n"%('\t'.join(data)))
                        
                        # Write new alleles
                        fout.write("%s\n"%('\t').join(new_allele))
                    
                    # Otherwise, copy alleles
                    else:
                        fout.write(line)
                        
    print("Finished")


if __name__ == "__main__":
    opts = get_opts()
    in_allele = opts.input
    min_num = opts.min_num
    out_allele = opts.output
    pull_down_paralog(in_allele, min_num, out_allele)
