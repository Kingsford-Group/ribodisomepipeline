#!/usr/bin/env python
import numpy as np
import scipy.stats
from footprint_hist_parser import parse_rlen_hist, get_cds_range
from ribomap_result_parser import parse_estimated_profile
from ribofit_utils import validate_profile, codonp_from_basep, sum_frames

def generate_cds_coverage_profile(tlist, rlen_min, rlen_max, cds_range):
    p = {}
    for rid in tlist:
        tid = tlist[rid]['tid']
        start, stop = cds_range[tid]
        profile = np.zeros(stop-start)
        for rlen, pos_cnt_list in tlist[rid]['prof'].iteritems():
            if rlen < rlen_min or rlen > rlen_max: 
                continue
            else:
                for pos, cnt in pos_cnt_list:
                    for i in xrange(pos-start, pos-start+rlen):
                        if i<0 or i>= stop-start: continue
                        profile[i] += cnt
        if np.any(profile!=0):
            p[tid] = profile
    return p

def get_total_disome_cnt(tlist, rlen_min, rlen_max, cds_range):
    tcnt = {}
    for rid in tlist:
        tid = tlist[rid]['tid']
        tot_cnt = 0
        start, stop = cds_range[tid]
        for rlen, pos_cnt_list in tlist[rid]['prof'].iteritems():
            if rlen < rlen_min or rlen > rlen_max:
                continue
            else:
                for pos, cnt in pos_cnt_list:
                    if pos < start or pos > stop : continue
                    tot_cnt += cnt
        if tot_cnt!=0:
            tcnt[tid] = tot_cnt
    return tcnt

def get_monocnt_from_ribomap_codon_prof(p):
    tcnt = {}
    for rid in p:
        tid = p[rid]['tid']
        tcnt[tid] = sum(p[rid]['rprofile'])
    return tcnt

def high_coverage_transcript_from_ribomap_codon_prof(p):
    return [ p[rid]['tid'] for rid in p if validate_profile(p[rid]['rprofile']) ]

def get_collision_rate(cds_range, doublet_fn, singlet_fn):
    doublets = parse_rlen_hist(doublet_fn)
     # singlets (dic of ribomap outputs: rid: tid, rprofile, nprofile, mprofile)
    sp_base = parse_estimated_profile(singlet_fn)
    sp_codon = codonp_from_basep(sp_base, cds_range, sum_frames)
    # high coverage transcript
    tid_list = high_coverage_transcript_from_ribomap_codon_prof(sp_codon)
    print "# high coverage transcript:", len(tid_list)
    # total ribo counts
    dcnt = get_total_disome_cnt(doublets, 58, 61, cds_range)
    print "# transcripts with doublets", len(dcnt)
    scnt = get_monocnt_from_ribomap_codon_prof(sp_codon)
    # collision rate
    cr = {}
    for tid in tid_list:
        if tid in dcnt:
            ratio = float(dcnt[tid])/(dcnt[tid] + scnt[tid])
            cr[tid] = ratio
    print "overlapped hc singlets and doublets: {0} ({1:.0%})".format(len(cr), len(cr)/float(len(tid_list)))
    print "average collision rate: {0:.2e}".format(np.mean(cr.values()))
    return cr

def correlate_two_sets(s1, s2):
    tid_list = set(s1.keys())&set(s2.keys())
    s1_common = [ s1[tid] for tid in tid_list ]
    s2_common = [ s2[tid] for tid in tid_list ]
    pr = scipy.stats.pearsonr(s1_common, s2_common)[0]
    return pr

if __name__ == "__main__":
    cds_txt = "../ref/cds_range.txt"
    nchx_dfn = "../rlen_hist/prelim/Lib-5-5-15_3_CACGAT_R1_nodup.hist"
    nchx_sfn = "../ribomap/Lib-5-5-15_2_AGTTCC_R1_nodup.base"
    nchx_nbc_dfn = "../rlen_hist/prelim/Lib-5-5-15_3_CACGAT_R1_nonempty.hist"
    nchx_nbc_sfn = "../ribomap/Lib-5-5-15_2_AGTTCC_R1_nonempty.base"
    chx_dfn = "../rlen_hist/prelim/Lib-5-5-15_7_ATCACG_R1_nodup.hist"
    chx_sfn = "../ribomap/Lib-5-5-15_6_TCCCGA_R1_nodup.base"
    chx_nbc_dfn = "../rlen_hist/prelim/Lib-5-5-15_7_ATCACG_R1_nonempty.hist"
    chx_nbc_sfn = "../ribomap/Lib-5-5-15_6_TCCCGA_R1_nonempty.base"
    wt_dfn = "../rlen_hist/Guydosh14/WT_disome.hist"
    wt_sfn = "../ribomap/
    cds_range = get_cds_range(cds_txt)
    cr_nchx = get_collision_rate(cds_range, nchx_dfn, nchx_sfn)
    cr_nchx_nbc = get_collision_rate(cds_range, nchx_nbc_dfn, nchx_nbc_sfn)    
    pr = correlate_two_sets(cr_nchx, cr_nchx_nbc)
    print "no Chx: correlation of barcode collapse before after {0}".format(pr)
    cr_chx = get_collision_rate(cds_range, chx_dfn, chx_sfn)
    cr_chx_nbc = get_collision_rate(cds_range, chx_nbc_dfn, chx_nbc_sfn)    
    pr = correlate_two_sets(cr_chx, cr_chx_nbc)
    print "Chx: correlation of barcode collapse before after {0}".format(pr)
    pr = correlate_two_sets(cr_chx_nbc, cr_nchx_nbc)
    print "before barcode collapse: correlate chx vs nchx {0}".format(pr)
    pr = correlate_two_sets(cr_chx, cr_nchx)
    print "after barcode collapse: correlate chx vs nchx {0}".format(pr)
