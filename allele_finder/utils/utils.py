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

        retain_db = {}
        dup_db = {}
        for _ in seq_db:
            if len(seq_db[_]) != 1:
                seq_ids = list(seq_db[_])
                random.shuffle(seq_ids)
                retain_id = seq_ids[0]
                for gid in seq_db[_]:
                    retain_db[gid] = retain_id
                dup_db[retain_id] = []
                for gid in seq_db[_]:
                    if gid != retain_id:
                        dup_db[retain_id].append(gid)

        return retain_db, dup_db

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

    @staticmethod
    def get_seq_length(in_fasta):
        seq_len_db = {}
        with open(in_fasta, 'r') as fin:
            for line in fin:
                if line[0] == '>':
                    gid = line.strip().split()[0][1:]
                    seq_len_db[gid] = 0
                else:
                    seq_len_db[gid] += len(line.strip())

        return seq_len_db


class BlastUtils:
    @staticmethod
    def allele_blast(in_blast, cds_len_db, threshold):
        allele_list = []
        used_genes = {}
        with open(in_blast, 'r') as fin:
            for line in fin:
                data = line.strip().split()
                qry = data[0]
                target = data[1]
                qry_len = cds_len_db[qry]
                target_len = cds_len_db[target]
                match = float(data[2]) / 100. * float(data[3])
                cur_iden = match * 2. / (qry_len + target_len)
                if cur_iden < threshold:
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
    def adjust_allele_table(in_allele, ref_gff3, ref_cds_len_db, hap_gff3, hap_cds_len_db,
                            ref_blast, hap_blast, blast_threshold, tandem, allele_num, out_allele,
                            is_mono, log_file):
        flog = open(log_file, 'w')
        flog.write("Loading ref blast\n")
        ref_blast_db = {}
        with open(ref_blast, 'r') as fin:
            for line in fin:
                data = line.strip().split()
                hap_gn = data[0]
                ref_gn = data[1]
                hap_len = hap_cds_len_db[hap_gn]
                ref_len = ref_cds_len_db[ref_gn]
                match = float(data[2]) / 100. * float(data[3])
                cur_iden = match * 2. / (hap_len + ref_len)
                if cur_iden < blast_threshold:
                    continue
                if hap_gn not in ref_blast_db or cur_iden > ref_blast_db[hap_gn][1]:
                    ref_blast_db[hap_gn] = [ref_gn, cur_iden]

        flog.write("Loading hap blast\n")
        hap_blast_db = {}
        with open(hap_blast, 'r') as fin:
            for line in fin:
                data = line.strip().split()
                gn1 = data[0]
                gn2 = data[1]
                glen1 = hap_cds_len_db[gn1]
                glen2 = hap_cds_len_db[gn2]
                match = float(data[2]) / 100. * float(data[3])
                cur_iden = match * 2. / (glen1 + glen2)
                if gn1 != gn2 and cur_iden >= blast_threshold:
                    if gn1 not in hap_blast_db or cur_iden > hap_blast_db[gn1][1]:
                        hap_blast_db[gn1] = [gn2, cur_iden]

        flog.write("Loading ref gff3\n")
        ref_db = {"NA": ["NA", "NA"]}
        with open(ref_gff3, 'r') as fin:
            for line in fin:
                if line.strip() == '' or line[0] == '#':
                    continue
                data = line.strip().split()
                if data[2] == 'gene':
                    gid = GFF3Utils.get_gene_id(data[8])
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
                gid = GFF3Utils.get_gene_id(data[8])
                hap_db[gid] = [chrn, int(data[3])]

        flog.write("Loading tandem\n")
        tandem_set = set()
        with open(tandem, 'r') as fin:
            for line in fin:
                data = line.strip().split(',')
                chrn_set = set()
                # tandem genes should in one chromosome
                for gid in data:
                    chrn_set.add(hap_db[gid][0])
                if len(chrn_set) == 1:
                    for gid in data:
                        tandem_set.add(gid)

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
                    if hidx < 0 or hidx >= allele_num:
                        continue
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
            new_info_db = {}
            for ref_gn in info_db:
                new_info_db[ref_gn] = ["" for _ in range(0, len(info_db[ref_gn]))]
                null_cnt = 0
                for i in range(0, len(info_db[ref_gn])):
                    if not info_db[ref_gn][i]:
                        new_info_db[ref_gn][i] = "NA"
                        null_cnt += 1
                    else:
                        tmp_list = sorted(info_db[ref_gn][i])
                        tmp_ids = [tmp_list[0][-1]]
                        if len(tmp_list) > 1:
                            for j in range(1, len(tmp_list)):
                                pos, gid = tmp_list[j]
                                if gid in tandem_set:
                                    tmp_ids.append(gid + "-T")
                                else:
                                    tmp_ids.append(gid + "-P")
                        new_info_db[ref_gn][i] = ','.join(tmp_ids)

                if null_cnt == allele_num:
                    continue
                if not is_mono or ref_gn.startswith("NA"):
                    ori_gn = ref_gn.split(":::")[0]
                else:
                    ori_gn = ref_gn
                tmp_list = copy.deepcopy(ref_db[ori_gn])
                tmp_list.append(ori_gn)
                tmp_list.extend(new_info_db[ref_gn])
                full_allele.append(tmp_list)

        flog.write("Writing allele table\n")
        with open(out_allele, 'w') as fout:
            allele_header = []
            for i in range(0, allele_num):
                allele_header.append("Allele %s" % (chr(i + 65)))
            fout.write("#CHR\tPOS\tRef gene\t%s\n" % ('\t'.join(allele_header)))
            for info in sorted(full_allele):
                fout.write("%s\n" % ('\t'.join(map(str, info))))

        flog.write("Finished\n")
        flog.close()


class GmapUtils:
    @staticmethod
    def read_gff3(in_gff3):
        drop_pre_set = {'ctg', 'utg', 'tig', 'sca', 'sup'}
        gff3_db = {}
        gene_order = {}
        gene_cnt = {}
        with open(in_gff3, 'r') as fin:
            for line in fin:
                if line.strip() == '' or line[0] == '#':
                    continue
                data = line.strip().split()
                if data[0][:3].lower() in drop_pre_set:
                    continue
                chrn = re.findall(r'([A-Za-z]+\d+).*', data[0])
                if chrn:
                    chrn = chrn[0]
                else:
                    continue
                feature_type = data[2]
                if feature_type != 'gene':
                    continue
                sp = int(data[3])
                ep = int(data[4])
                gene = GFF3Utils.get_gene_id(data[8])
                if chrn not in gff3_db:
                    gff3_db[chrn] = []
                if chrn not in gene_cnt:
                    gene_cnt[chrn] = 0
                gene_cnt[chrn] += 1
                gene_order[gene] = [chrn, gene_cnt[chrn]]
                gff3_db[chrn].append([sp, ep, gene])
        return gff3_db, gene_order

    @staticmethod
    def merge_allele(region_list, overlap_ratio):
        alleles = []
        last_ep = 0
        last_len = 0
        for sp, ep, gene in sorted(region_list):
            cur_len = ep - sp + 1
            if last_ep == 0:
                alleles.append([])
                alleles[-1].append(gene)
            else:
                if sp > last_ep or last_ep - sp + 1 < ((cur_len + last_len) / 2.) * overlap_ratio:
                    alleles.append([])
                alleles[-1].append(gene)
            if ep > last_ep:
                last_ep = ep
                last_len = cur_len
        return alleles

    @staticmethod
    def allele_gmap(gff3_db, overlap_ratio, threads):
        pool = multiprocessing.Pool(processes=threads)
        res = {}
        for chrn in gff3_db:
            res[chrn] = pool.apply_async(GmapUtils.merge_allele, (gff3_db[chrn], overlap_ratio,))
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
    def bin_search(qry_pos, pos_list):
        s = 0
        e = len(pos_list) - 1
        target_idx = -1
        while s <= e:
            mid = int((s + e) / 2)
            if pos_list[mid][0] < qry_pos:
                s = mid + 1
            elif pos_list[mid][0] > qry_pos:
                e = mid - 1
            else:
                target_idx = mid
                break
        return target_idx, e

    @staticmethod
    def check_ovlp(sp, ep, pos_list, thres):
        lbound, candidate_lbound = TEUtils.bin_search(sp, pos_list)
        if lbound == -1:
            if candidate_lbound == -1:
                candidate_lbound = 0
            if pos_list[candidate_lbound][1] >= sp:
                lbound = candidate_lbound
            else:
                lbound = candidate_lbound + 1
        rbound, candidate_rbound = TEUtils.bin_search(ep, pos_list)
        if rbound == -1:
            if candidate_rbound == -1:
                candidate_rbound = 0
            rbound = candidate_rbound

        ovlp_len = 0
        if lbound <= rbound:
            tmp_regions = []
            for rsp, rep in pos_list[lbound: rbound + 1]:
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
    def filter_with_TE(in_allele, hap_gff3, in_te, te_thres, te_filter_only_paralog, out_allele, log_file):
        flog = open(log_file, 'w')
        flog.write("Loading TEs\n")
        TE_db = {}
        with open(in_te, 'r') as fin:
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
        ovlp_ids = set()
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
                    if TEUtils.check_ovlp(sp, ep, TE_db[chrn], te_thres):
                        ovlp_ids.add(gid)

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
                            if not te_filter_only_paralog:
                                for gid in ids:
                                    if '-' in gid:
                                        rid = gid.split('-')[0]
                                    else:
                                        rid = gid
                                    if rid in ovlp_ids:
                                        remove_ids.append(rid)
                                        continue
                                    tmp_ids.append(gid)
                                if not tmp_ids:
                                    tmp_ids.append('NA')
                                else:
                                    normal_ids = []
                                    tandem_ids = []
                                    paralog_ids = []
                                    for gid in tmp_ids:
                                        if '-' not in gid:
                                            normal_ids.append(gid)
                                        else:
                                            if gid[-1] == 'P':
                                                paralog_ids.append(gid)
                                            else:
                                                tandem_ids.append(gid)
                                    if not normal_ids:
                                        if paralog_ids:
                                            normal_ids = [paralog_ids[0].split('-')[0]]
                                            paralog_ids = paralog_ids[1:]
                                        elif tandem_ids:
                                            normal_ids = [tandem_ids[0].split('-')[0]]
                                            tandem_ids = tandem_ids[1:]
                                        else:
                                            normal_ids = []
                                    tmp_ids = normal_ids
                                    tmp_ids.extend(tandem_ids)
                                    tmp_ids.extend(paralog_ids)
                            else:
                                for gid in ids:
                                    if '-' in gid and gid[-1] == 'P':
                                        rid = gid.split('-')[0]
                                        if rid in ovlp_ids:
                                            remove_ids.append(rid)
                                            continue
                                    tmp_ids.append(gid)

                            data[i] = ','.join(tmp_ids)
                        # drop all NA line
                        if len(set(data[3:])) == 1:
                            continue
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
        last_line = ""
        with open(in_gff3, 'r') as fin:
            with open(out_gff3, 'w') as fout:
                for line in fin:
                    if line.strip() == '' or line[0] == '#':
                        if line == last_line:
                            continue
                        fout.write(line)
                        last_line = line
                    else:
                        data = line.strip().split()
                        if data[2] == 'gene':
                            if "Name" in data[8]:
                                gid = re.findall(r'Name=(.*)', data[8])[0].split(';')[0]
                            else:
                                gid = re.findall(r'ID=(.*)', data[8])[0].split(';')[0]

                        if gid not in drop_id_set:
                            fout.write(line)
                            last_line = line

    @staticmethod
    def get_gene_id(gff_id_col):
        if "Name" in gff_id_col:
            gid = re.findall(r'Name=(.*)', gff_id_col)[0].split(';')[0]
        else:
            gid = re.findall(r'ID=(.*)', gff_id_col)[0].split(';')[0]
        return gid
