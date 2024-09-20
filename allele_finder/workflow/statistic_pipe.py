import re
from allele_finder.utils.message import Message


def stat_formatted_allele(allele_file, hap_gff3, out_pre):
    Message.info("Loading hap gff3")
    hap_db = {}
    with open(hap_gff3, 'r') as fin:
        for line in fin:
            if line.strip() == '' or line[0] == '#':
                continue
            data = line.strip().split()
            if 'ctg' in data[0] or 'tig' in data[0] or 'utg' in data[0]:
                continue
            chrn = "Chr%02d" % (int(re.findall(r'[A-Za-z](\d+)', data[0])[0]))
            if "Name" in data[8]:
                gid = re.findall(r'Name=(.*)', data[8])[0].split(';')[0]
            else:
                gid = re.findall(r'ID=(.*)', data[8])[0].split(';')[0]
            hap_db[gid] = chrn

    Message.info("Loading allele table")
    allele_summary = {}
    gene_cnt_summary = {}
    with open(allele_file, 'r') as fin:
        for line in fin:
            data = line.strip().split('\t')
            if line[0] == '#':
                allele_cnt = len(data) - 3
                continue

            tand_cnt = 0
            disp_cnt = 0
            allele_cnt_list = [0 for i in range(0, allele_cnt)]
            gene_cnt = 0
            for i in range(3, len(data)):
                if data[i] == 'NA':
                    continue
                gene_cnt += 1
                ids = data[i].split(',')
                chrn_cnt_db = {}
                allele_cnt_list[i - 3] += 1
                for gid in ids:
                    if '-' in gid and gid.split('-')[-1] in {"T", "P"}:
                        tmp = gid.split('-')
                        repeat_type = tmp[-1]
                        gid = '-'.join(tmp[:-1])
                        # id, type = id.split('-')
                        if repeat_type == 'T':
                            tand_cnt += 1
                        else:
                            disp_cnt += 1
                    if gid in hap_db:
                        chrn = hap_db[gid]
                        if chrn not in chrn_cnt_db:
                            chrn_cnt_db[chrn] = 0
                        chrn_cnt_db[chrn] += 1
            max_chrn = ""
            max_cnt = 0
            for chrn in chrn_cnt_db:
                if chrn_cnt_db[chrn] > max_cnt:
                    max_cnt = chrn_cnt_db[chrn]
                    max_chrn = chrn
            if max_chrn not in allele_summary:
                # [[allele count list], tandem count, dispersely count, line count]
                allele_summary[max_chrn] = [[0 for i in range(0, allele_cnt)], 0, 0, 0]

            allele_summary[max_chrn][-1] += 1
            for i in range(0, allele_cnt):
                allele_summary[max_chrn][0][i] += allele_cnt_list[i]

            allele_summary[max_chrn][1] += tand_cnt
            allele_summary[max_chrn][2] += disp_cnt

            if max_chrn not in gene_cnt_summary:
                gene_cnt_summary[max_chrn] = [0 for i in range(0, allele_cnt)]
            gene_cnt_summary[max_chrn][gene_cnt - 1] += 1

    Message.info("Writing summary")
    with open(out_pre + '.allele.stat', 'w') as fout:
        allele_list = []
        for i in range(0, allele_cnt):
            allele_list.append("Allele %s" % chr(65 + i))
        fout.write("#CHR,Line count,%s,Tandem,Dispersely\n" % (','.join(allele_list)))
        sum_info = [0 for i in range(0, allele_cnt + 3)]
        for chrn in sorted(allele_summary):
            info = [chrn, allele_summary[chrn][3]]
            info.extend(allele_summary[chrn][0])
            info.append(allele_summary[chrn][1])
            info.append(allele_summary[chrn][2])
            sum_info[0] += allele_summary[chrn][3]
            for i in range(0, len(allele_summary[chrn][0])):
                sum_info[i + 1] += allele_summary[chrn][0][i]
            sum_info[allele_cnt + 1] += allele_summary[chrn][1]
            sum_info[allele_cnt + 2] += allele_summary[chrn][2]
            fout.write("%s\n" % (','.join(map(str, info))))
        fout.write("Total,%s\n" % (','.join(map(str, sum_info))))

    with open(out_pre + '.genes.stat', 'w') as fout:
        gene_cnt_list = []
        sum_info = [0 for i in range(0, allele_cnt + 3)]
        for i in range(allele_cnt, 0, -1):
            gene_cnt_list.append('No. with %d' % i)
        fout.write("#CHR,Total count of allele genes,%s,Tandem,Dispersely\n" % (','.join(gene_cnt_list)))
        for chrn in sorted(gene_cnt_summary):
            info = [chrn, sum(gene_cnt_summary[chrn])]
            for i in range(allele_cnt, 0, -1):
                info.append(gene_cnt_summary[chrn][i - 1])
            info.append(allele_summary[chrn][1])
            info.append(allele_summary[chrn][2])
            for i in range(0, len(sum_info)):
                sum_info[i] += info[i + 1]
            fout.write("%s\n" % (','.join(map(str, info))))
        fout.write("Total,%s\n" % (','.join(map(str, sum_info))))

    Message.info("Finished")


def main(args):
    allele_file = args.input
    hap_gff3 = args.gff3
    out_pre = args.output
    stat_formatted_allele(allele_file, hap_gff3, out_pre)
