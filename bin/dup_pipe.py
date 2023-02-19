#!/usr/bin/env python
import copy

def read_dup_class(gene_type_file):
    """
    Duplicate 1	Location	Duplicate 2	Location	E-value
    AT1G01010.1	Ath-Chr1:3631	AT4G01550.1	Ath-Chr4:673862	5e-52
    AT1G01020.1	Ath-Chr1:5928	AT4G01510.1	Ath-Chr4:642733	5e-74
    AT1G01030.1	Ath-Chr1:11649	AT1G13260.1	Ath-Chr1:4542168	3e-45
    """
    with open(gene_type_file,'r') as fin:
        dup_relation = {}
        dup_lst = []
        fin.readline()
        for line in fin:
            tmpLst = line.split('\t')
            g1, g2 = tmpLst[0], tmpLst[2]
            if g1 not in dup_relation:
                dup_relation[g1] = []
            dup_relation[g1].append(g2)
            dup_lst.append(g2)
    return dup_relation, dup_lst

def merge_dup_class(outputdir):
    tandemDic, tandemLst = read_dup_class(outputdir + "xyz.tandem.pairs")
    proximalDic, proximalLst = read_dup_class(outputdir + "xyz.proximal.pairs")
    transREpDic, transREpLst = read_dup_class(outputdir + "xyz.transposed.pairs")
    dispersedDic, dispersedLst = read_dup_class(outputdir + "xyz.dispersed.pairs")
    paralogDic = {**proximalDic, **transREpDic, **dispersedDic}
    paralogLst = [*proximalLst, *transREpLst, *dispersedLst]
    return paralogDic, paralogLst, tandemDic, tandemLst

def remove_repeat_allele(base_allele, pLst, tLst):
    new_allele_lst = []
    for allele_list in base_allele:
        tmpLst = []
        for allele in allele_list:
            if allele in pLst or allele in tLst:
                continue
            else:
                tmpLst.append(allele)
        if tmpLst != []:
            new_allele_lst.append(tmpLst)
        else:
            pass
    return new_allele_lst

                
def add_dup_info(full_allele, info_db, hapdb, singleton, allele2table, addition_allele, pDic, tDic):
    for tmpItemIndex in range(len(full_allele)):
        tmpItem = full_allele[tmpItemIndex]
        ref_gn = tmpItem[2]
        if ref_gn in allele2table:
            refInfo = tmpItem[:3]
            alleleInfo = info_db[ref_gn]
            refInfo.extend(alleleInfo)
            full_allele[tmpItemIndex] = refInfo
    for sg in singleton:
        ori_pos = copy.deepcopy(hapdb[sg])
        ori_pos.append(sg)
        ori_pos.extend(info_db[sg])
        full_allele.append(ori_pos)
    for aa in addition_allele:
        ori_pos_aa = copy.deepcopy(hapdb[aa])
        ori_pos_aa.append(aa)
        ori_pos_aa.extend(info_db[aa])
        full_allele.append(ori_pos_aa)
    # reform allele table
    for tmpItemIndex in range(len(full_allele)):
        tmpItem = full_allele[tmpItemIndex]
        tmpAlleleInfo = tmpItem[3:]
        tmpAlleleInfo = [[tmpAlleleInfo[i]] for i in range(len(tmpAlleleInfo))]
        for aI in range(len(tmpAlleleInfo)):
            tmpout = copy.deepcopy(tmpAlleleInfo[aI])
            for g in tmpAlleleInfo[aI]:
                if g in pDic:
                    for rg in pDic[g]:
                        tmpParalog = pDic[rg] + '-P'
                        tmpout.append(tmpParalog)
                if g in tDic:
                    for rg in tDic[g]:
                        tmpTandem = tDic[g] + '-T'
                        tmpout.append(tmpTandem)
            tmpout = ','.join(tmpout)
            tmpAlleleInfo[aI] = tmpout
        tmpRefInfo = tmpItem[:3]
        tmpRefInfo.extend(tmpAlleleInfo)
        full_allele[tmpItemIndex] = tmpRefInfo
        """
        full_allele -> [[ref_Chr, ref_pos, ref_gn, hapA, hapB, hapC, ..., hapE], [...], ...]
        pDic -> {gene1: gene1-P}
        pLst -> [gene1-P]
        """
    return full_allele

def read_repeat_file(rf):
    rDic = {}
    with open(rf, 'r') as fin:
        for line in fin:
            tmp1 = line.rstrip().split('\t')
            if tmp1[0] not in rDic:
                rDic[tmp1[0]] = []
            tmp2 = tmp1[1].split(',')
            for t in tmp2:
                rDic[tmp1[0]].append(t)
    return rDic

def write_repeat_file(rDic, rpath):
    with open(rpath, 'w') as fout:
        for r in rDic:
            fout.write("{}\t{}\n".format(r, ','.join(rDic[r])))
