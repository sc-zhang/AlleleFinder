import re
import multiprocessing
import copy
import random


class FastaUtils:
    @staticmethod
    def get_full_same_seqs(fasta_file):
        seq_db = {}
        with open(fasta_file, 'r') as fin:
            gid = ""
            seq = []
            for line in fin:
                if line.strip() == "":
                    continue
                if line[0] == '>':
                    if seq:
                        full_seq = ''.join(seq)
                        if full_seq not in seq_db:
                            seq_db[full_seq] = set()
                        seq_db[full_seq].add(gid)
                    gid = line.strip()[1:]
                    seq = []
                else:
                    seq.append(line.strip().upper())
            if seq:
                full_seq = ''.join(seq)
                if full_seq not in seq_db:
                    seq_db[full_seq] = set()
                seq_db[full_seq].add(gid)

        dup_db = {}
        for _ in seq_db:
            if len(seq_db[_]) != 1:
                seq_ids = list(seq_db[_])
                random.shuffle(seq_ids)
                retain_id = seq_ids[0]
                for gid in seq_db[_]:
                    dup_db[gid] = retain_id

        return dup_db

    @staticmethod
    def filter_fasta(in_fasta, out_fasta, drop_id_set):
        fa_db = {}
        with open(in_fasta, 'r') as fin:
            for line in fin:
                if line.strip() == '':
                    continue
                if line[0] == '>':
                    gid = line.strip()[1:]
                    fa_db[gid] = []
                else:
                    fa_db[gid].append(line.strip().upper())
        for gid in fa_db:
            fa_db[gid] = ''.join(fa_db[gid])

        with open(out_fasta, 'w') as fout:
            for gid in sorted(fa_db):
                if gid in drop_id_set:
                    continue
                fout.write(">%s\n%s\n" % (gid, fa_db[gid]))

    @staticmethod
    def check_fasta_type(fasta_file):
        with open(fasta_file, 'r') as fin:
            for line in fin:
                if line[0] == ">":
                    continue
                if "M" in line.upper():
                    return "prot"
                else:
                    return "nucl"

    @staticmethod
    def get_seq_ids(fasta_file):
        gene_id_set = set()
        with open(fasta_file, 'r') as fin:
            for line in fin:
                if not line.strip():
                    continue
                if line[0] == '>':
                    gene_id_set.add(line.strip().split()[0][1:])
        return gene_id_set

    @staticmethod
    def split_fasta_with_allele(in_fa, in_allele, out_pre):
        fa_db = {}
        with open(in_fa, 'r') as fin:
            for line in fin:
                if line[0] == '>':
                    gid = line.strip().split()[0][1:]
                    fa_db[gid] = []
                else:
                    fa_db[gid].append(line.strip())

        for gid in fa_db:
            fa_db[gid] = ''.join(fa_db[gid])

        single_list = []
        multi_list = []
        with open(in_allele, 'r') as fin:
            for line in fin:
                data = line.strip().split(',')
                cnt = 0
                for gene in data:
                    if gene.strip() != '':
                        cnt += 1
                if cnt == 1:
                    single_list.append(data[0])
                else:
                    for gene in data:
                        if gene.strip() != '':
                            multi_list.append(gene)

        with open(out_pre + '_single.fa', 'w') as fout:
            for gid in single_list:
                fout.write(">%s\n%s\n" % (gid, fa_db[gid]))

        with open(out_pre + '_multi.fa', 'w') as fout:
            for gid in multi_list:
                fout.write(">%s\n%s\n" % (gid, fa_db[gid]))


class BlastUtils:
    @staticmethod
    def allele_blast(in_blast, iden):
        allele_list = []
        used_genes = {}
        with open(in_blast, 'r') as fin:
            for line in fin:
                data = line.strip().split()
                qry = data[0]
                target = data[1]
                if float(data[2]) < iden:
                    continue
                if qry not in used_genes:
                    used_genes[qry] = 1
                    allele_list.append([qry, target])
        return allele_list


class AlleleUtils:

    @staticmethod
    def merge_allele(allele_list):
        pos_db = {}
        new_allele_list = []
        for i in range(0, len(allele_list)):
            merge_list = []
            scrappy_allele = []
            for allele in allele_list[i]:
                if allele in pos_db:
                    merge_list.append(allele)
                else:
                    scrappy_allele.append(allele)
            if not merge_list:
                new_allele_list.append(allele_list[i])
                for allele in allele_list[i]:
                    pos_db[allele] = len(new_allele_list) - 1
            else:
                merge_pos = pos_db[merge_list[0]]
                if scrappy_allele:
                    for allele in scrappy_allele:
                        pos_db[allele] = merge_pos
                    new_allele_list[merge_pos].extend(scrappy_allele)
                for allele in merge_list[1:]:
                    allele_pos = pos_db[allele]
                    if allele_pos == merge_pos:
                        continue
                    new_allele_list[merge_pos].extend(new_allele_list[allele_pos])
                    for other_allele in new_allele_list[allele_pos]:
                        pos_db[other_allele] = merge_pos
                    new_allele_list[allele_pos] = []
        allele_list = []
        for i in range(0, len(new_allele_list)):
            if new_allele_list[i]:
                allele_list.append(list(set(new_allele_list[i])))

        return allele_list

    @staticmethod
    def split_allele(allele_list, gene_db):
        new_allele_list = []
        for i in range(0, len(allele_list)):
            tmp_list = allele_list[i]
            tmp_db = {}
            for gene in tmp_list:
                if gene not in gene_db:
                    continue
                chrn = gene_db[gene][0]
                if chrn not in tmp_db:
                    tmp_db[chrn] = []
                tmp_db[chrn].append(gene)
            for chrn in tmp_db:
                new_allele_list.append(tmp_db[chrn])
        return new_allele_list

    @staticmethod
    def get_allele_with_mcscanx(html_file):
        allele_list = []
        with open(html_file, 'r') as fin:
            for line in fin:
                if line.startswith("<tr"):
                    vals = re.findall('<td.*?>(.*?)</td>', line)
                    tmp_list = []
                    if len(vals) < 2:
                        continue
                    for val in vals[1:]:
                        if '&' not in val and 'ctg' not in val:
                            tmp_list.append(val)
                    if tmp_list:
                        allele_list.append(tmp_list)
        return allele_list

    @staticmethod
    def adjust_allele_table(in_allele, ref_gff3, hap_gff3, ref_blast, hap_blast, iden, tandem, allele_num, out_allele,
                            is_mono, log_file):
        flog = open(log_file, 'w')
        flog.write("Loading ref blast\n")
        ref_blast_db = {}
        with open(ref_blast, 'r') as fin:
            for line in fin:
                data = line.strip().split()
                hap_gn = data[0]
                ref_gn = data[1]
                bi = float(data[2])
                if bi < iden:
                    continue
                bs = data[-1]
                if hap_gn not in ref_blast_db:
                    ref_blast_db[hap_gn] = [ref_gn, bs]
                if bs > ref_blast_db[hap_gn][-1]:
                    ref_blast_db[hap_gn] = [ref_gn, bs]

        flog.write("Loading hap blast\n")
        hap_blast_db = {}
        with open(hap_blast, 'r') as fin:
            for line in fin:
                data = line.strip().split()
                gn1 = data[0]
                gn2 = data[1]
                bi = float(data[2])
                if gn1 != gn2 and bi >= iden:
                    if gn1 not in hap_blast_db:
                        hap_blast_db[gn1] = [gn2, bi]
                    elif bi > hap_blast_db[gn1][1]:
                        hap_blast_db[gn1] = [gn2, bi]

        flog.write("Loading ref gff3\n")
        ref_db = {"NA": ["NA", "NA"]}
        with open(ref_gff3, 'r') as fin:
            for line in fin:
                if line.strip() == '' or line[0] == '#':
                    continue
                data = line.strip().split()
                if data[2] == 'gene':
                    if "Name" in data[8]:
                        gid = re.findall(r'Name=(.*)', data[8])[0].split(';')[0]
                    else:
                        gid = re.findall(r'ID=(.*)', data[8])[0].split(';')[0]
                    if gid not in ref_db:
                        ref_db[gid] = [data[0], int(data[3])]

        flog.write("Loading hap gff3\n")
        hap_db = {}
        with open(hap_gff3, 'r') as fin:
            for line in fin:
                if line.strip() == '' or line[0] == '#':
                    continue
                data = line.strip().split()
                if data[2] != 'gene':
                    continue
                chrn = data[0]
                if "Name" in data[8]:
                    gid = re.findall(r'Name=(.*)', data[8])[0].split(';')[0]
                else:
                    gid = re.findall(r'ID=(.*)', data[8])[0].split(';')[0]
                hap_db[gid] = [chrn, int(data[3])]

        flog.write("Loading tandem\n")
        tandem_db = {}
        with open(tandem, 'r') as fin:
            for line in fin:
                data = line.strip().split(',')
                for gid in data:
                    tandem_db[gid] = 1

        flog.write("Formatting allele table\n")
        full_allele = []
        info_db = {}
        line_cnt = 0
        with open(in_allele, 'r') as fin:
            for line in fin:
                data = line.strip().split(',')
                line_cnt += 1
                for gid in data:
                    if gid not in hap_db:
                        continue
                    hchrn, hpos = hap_db[gid]
                    hidx = ord(hchrn[-1]) - 65
                    if gid not in ref_blast_db:
                        if gid not in hap_blast_db:
                            ref_id = "NA"
                        else:
                            ref_id = hap_blast_db[gid][0]
                    else:
                        ref_id = gid
                    if ref_id not in ref_blast_db or ref_id not in data:
                        ref_gn = "NA:::%d" % line_cnt
                    else:
                        ref_gn = ref_blast_db[ref_id][0]
                    if is_mono:
                        if ref_gn not in info_db:
                            info_db[ref_gn] = [[] for _ in range(0, allele_num)]
                        info_db[ref_gn][hidx].append([hpos, gid])
                    else:
                        tmp_gn = "%s:::%d" % (ref_gn, line_cnt)
                        if tmp_gn not in info_db:
                            info_db[tmp_gn] = [[] for _ in range(0, allele_num)]
                        info_db[tmp_gn][hidx].append([hpos, gid])

            for ref_gn in info_db:
                null_cnt = 0
                for i in range(0, len(info_db[ref_gn])):
                    if not info_db[ref_gn][i]:
                        info_db[ref_gn][i] = "NA"
                        null_cnt += 1
                    else:
                        tmp_list = sorted(info_db[ref_gn][i])
                        tmp_ids = [tmp_list[0][-1]]
                        if len(tmp_list) > 1:
                            for j in range(1, len(tmp_list)):
                                pos, gid = tmp_list[j]
                                if gid in tandem_db:
                                    tmp_ids.append(gid + "-T")
                                else:
                                    tmp_ids.append(gid + "-P")
                        info_db[ref_gn][i] = ','.join(tmp_ids)

                if null_cnt == allele_num:
                    continue
                if not is_mono or ref_gn.startswith("NA"):
                    ori_gn = ref_gn.split(":::")[0]
                else:
                    ori_gn = ref_gn
                tmp_list = copy.deepcopy(ref_db[ori_gn])
                tmp_list.append(ori_gn)
                tmp_list.extend(info_db[ref_gn])
                full_allele.append(tmp_list)

        flog.write("Writing allele table\n")
        with open(out_allele, 'w') as fout:
            allele_header = []
            for i in range(0, allele_num):
                allele_header.append("Allele %s" % (chr(i + 65)))
            fout.write("#CHR\tPOS\tref gene\t%s\n" % ('\t'.join(allele_header)))
            for info in sorted(full_allele):
                fout.write("%s\n" % ('\t'.join(map(str, info))))

        flog.write("Finished\n")
        flog.close()


class GmapUtils:
    @staticmethod
    def read_gff3(in_gff3):
        gff3_db = {}
        gene_order = {}
        gene_cnt = {}
        with open(in_gff3, 'r') as fin:
            for line in fin:
                if line.strip() == '' or line[0] == '#':
                    continue
                data = line.strip().split()
                if 'ctg' in data[0] or 'tig' in data[0] or 'utg' in data[0] or 'scaffold' in data[0] or 'super' in data[
                    0]:
                    continue
                chrn = re.findall(r'([A-Za-z]+\d+).*', data[0])[0]
                feature_type = data[2]
                if feature_type != 'gene':
                    continue
                sp = int(data[3])
                ep = int(data[4])
                if 'Name' in data[8]:
                    gene = re.findall(r'Name=(.*)', data[8])[0].split(';')[0]
                else:
                    gene = re.findall(r'ID=(.*)', data[8])[0].split(';')[0]
                if 'tig' in gene or 'ctg' in gene:
                    continue
                if chrn not in gff3_db:
                    gff3_db[chrn] = []
                if chrn not in gene_cnt:
                    gene_cnt[chrn] = 0
                gene_cnt[chrn] += 1
                gene_order[gene] = [chrn, gene_cnt[chrn]]
                gff3_db[chrn].append([sp, ep, gene])
        return gff3_db, gene_order

    @staticmethod
    def merge_allele(region_list):
        alleles = []
        last_ep = 0
        for sp, ep, gene in sorted(region_list):
            if last_ep == 0:
                alleles.append([])
                alleles[-1].append(gene)
                last_ep = ep
            else:
                if sp > last_ep:
                    alleles.append([])
                alleles[-1].append(gene)
                last_ep = ep
        return alleles

    @staticmethod
    def allele_gmap(gff3_db, threads):
        pool = multiprocessing.Pool(processes=threads)
        res = {}
        for chrn in gff3_db:
            res[chrn] = pool.apply_async(GmapUtils.merge_allele, (gff3_db[chrn],))
        pool.close()
        pool.join()
        allele_list = []
        for chrn in res:
            allele_list.extend(res[chrn].get())
        return allele_list


class TEUtils:
    @staticmethod
    def merge_regions(regions):
        tmp_regions = []
        last_ep = 0
        for sp, ep in sorted(regions):
            if not tmp_regions:
                tmp_regions.append(sp)
                last_ep = ep
            else:
                if sp > last_ep:
                    tmp_regions.append(last_ep)
                    tmp_regions.append(sp)
                    last_ep = ep
                else:
                    if ep > last_ep:
                        last_ep = ep
        tmp_regions.append(last_ep)
        new_regions = []
        for i in range(0, len(tmp_regions) - 1, 2):
            new_regions.append([tmp_regions[i], tmp_regions[i + 1]])
        return new_regions

    @staticmethod
    def check_ovlp(gid, sp, ep, pos_list, thres):
        s = 0
        e = len(pos_list) - 1
        l = -1
        while s <= e:
            mid = int((s + e) / 2)
            if pos_list[mid][0] < sp:
                s = mid + 1
            elif pos_list[mid][0] > sp:
                e = mid - 1
            else:
                l = mid
                break
        if l == -1:
            if e == -1:
                e = 0
            if pos_list[e][1] >= sp:
                l = e
            else:
                l = e + 1
        s = 0
        e = len(pos_list) - 1
        r = -1
        while s <= e:
            mid = int((s + e) / 2)
            if pos_list[mid][0] < ep:
                s = mid + 1
            elif pos_list[mid][0] > ep:
                e = mid - 1
            else:
                r = mid
                break
        if r == -1:
            if e == -1:
                e = 0
            r = e

        ovlp_len = 0
        if l <= r:
            tmp_regions = []
            for rsp, rep in pos_list[l: r + 1]:
                ovlp = min(ep, rep) - max(sp, rsp) + 1
                if ovlp <= 0:
                    continue
                tmp_regions.append([max(sp, rsp), min(ep, rep)])

            if len(tmp_regions) != 0:
                tmp_regions = TEUtils.merge_regions(tmp_regions)
                for msp, mep in tmp_regions:
                    ovlp_len += mep - msp + 1

        if ovlp_len * 1.0 / (ep - sp + 1) > thres:
            return True
        else:
            return False

    @staticmethod
    def filter_with_TE(in_allele, hap_gff3, in_TE, TE_thres, out_allele, log_file):
        flog = open(log_file, 'w')
        flog.write("Loading TEs\n")
        TE_db = {}
        with open(in_TE, 'r') as fin:
            for line in fin:
                if line.strip() == '' or line[0] == '#':
                    continue
                data = line.strip().split()
                chrn = data[0]
                sp = int(data[3])
                ep = int(data[4])

                if chrn not in TE_db:
                    TE_db[chrn] = []
                TE_db[chrn].append([sp, ep])
        for chrn in TE_db:
            TE_db[chrn] = TEUtils.merge_regions(sorted(TE_db[chrn]))

        flog.write("Getting overlap genes with TE\n")
        ovlp_ids = {}
        with open(hap_gff3, 'r') as fin:
            for line in fin:
                if line.strip() == '' or line[0] == '#':
                    continue
                data = line.strip().split()
                if data[2] != 'gene':
                    continue
                chrn = data[0]
                sp = int(data[3])
                ep = int(data[4])
                if 'Name' in data[8]:
                    gid = re.findall(r'Name=(.*)', data[8])[0].split(';')[0]
                else:
                    gid = re.findall(r'ID=(.*)', data[8])[0].split(';')[0]

                if chrn in TE_db:
                    if TEUtils.check_ovlp(gid, sp, ep, TE_db[chrn], TE_thres):
                        ovlp_ids[gid] = 1

        flog.write("Filtering allele table\n")
        remove_ids = []
        with open(in_allele, 'r') as fin:
            with open(out_allele, 'w') as fout:
                for line in fin:
                    if line[0] == '#':
                        fout.write(line)
                    else:
                        data = line.strip().split()
                        for i in range(3, len(data)):
                            if data[i] == 'NA':
                                continue
                            ids = data[i].split(',')
                            tmp_ids = []
                            for gid in ids:
                                if '-' in gid and gid[-1] == 'P':
                                    rid = gid.split('-')[0]
                                    if rid in ovlp_ids:
                                        remove_ids.append(rid)
                                        continue
                                tmp_ids.append(gid)
                            data[i] = ','.join(tmp_ids)
                        fout.write("%s\n" % ('\t'.join(data)))
        remove_fn = out_allele.split('.')
        remove_fn[-1] = 'removed.txt'
        remove_fn = '.'.join(remove_fn)
        with open(remove_fn, 'w') as fout:
            fout.write('\n'.join(sorted(remove_ids)))
        flog.write("Finished\n")
        flog.close()


class GFF3Utils:
    @staticmethod
    def filter_gff3(in_gff3, out_gff3, drop_id_set):
        with open(in_gff3, 'r') as fin:
            with open(out_gff3, 'w') as fout:
                is_write = False
                for line in fin:
                    if line.strip() == '' or line[0] == '#':
                        fout.write(line)
                    else:
                        data = line.strip().split()
                        if data[2] == 'gene':
                            if "Name" in data[8]:
                                gid = re.findall(r'Name=(.*)', data[8])[0].split(';')[0]
                            else:
                                gid = re.findall(r'ID=(.*)', data[8])[0].split(';')[0]
                            if gid in drop_id_set:
                                is_write = False
                            else:
                                is_write = True
                        if is_write:
                            fout.write(line)
