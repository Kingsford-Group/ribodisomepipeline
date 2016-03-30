#!/usr/bin/env python
import os
import numpy as np
import scipy.stats
from footprint_hist_parser import parse_rlen_hist, get_cds_range
from ribomap_result_parser import parse_estimated_profile
from ribofit_utils import validate_profile, codonp_from_basep, sum_frames
from file_names import *

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
    print os.path.basename(singlet_fn)
    doublets = parse_rlen_hist(doublet_fn)
     # singlets (dic of ribomap outputs: rid: tid, rprofile, nprofile, mprofile)
    sp_base = parse_estimated_profile(singlet_fn)
    sp_codon = codonp_from_basep(sp_base, cds_range, sum_frames)
    # total ribo counts
    dcnt = get_total_disome_cnt(doublets, dlen_min, dlen_max, cds_range)
    print "# transcripts with doublets", len(dcnt)
    scnt = get_monocnt_from_ribomap_codon_prof(sp_codon)
    # high coverage transcript
    # tid_list = high_coverage_transcript_from_ribomap_codon_prof(sp_codon)
    tid_list = list(set(dcnt.keys())&set(scnt.keys()))
    print "# high coverage transcript:", len(tid_list)
    # collision rate
    cr = {}
    for tid in tid_list:
        if tid in dcnt:
            ratio = float(dcnt[tid])/(dcnt[tid] + scnt[tid])
            cr[tid] = ratio
    print "high-covered transcripts with non-zero doublets: {0} ({1:.0%})".format(len(cr), len(cr)/float(len(tid_list)))
    print "average collision rate: {0:.2%}\n".format(np.mean(cr.values()))
    # build singlet sum and doublet sum as well
    shc = { tid:scnt[tid] for tid in cr }
    dhc = { tid:dcnt[tid] for tid in cr }
    return cr, shc, dhc

def correlate_two_sets(s1, s2):
    tid_list = set(s1.keys())&set(s2.keys())
    s1_common = [ s1[tid] for tid in tid_list ]
    s2_common = [ s2[tid] for tid in tid_list ]
    pr = scipy.stats.pearsonr(s1_common, s2_common)[0]
    return pr

if __name__ == "__main__":
    cds_range = get_cds_range(cds_txt)

    cr_nchx, scnt_nchx, dcnt_nchx = get_collision_rate(cds_range, nchx_dfn, nchx_sfn)
    # cr_nchx_nbc, scnt_nchx_nbc, dcnt_nchx_nbc = get_collision_rate(cds_range, nchx_nbc_dfn, nchx_nbc_sfn)    
    # pr = correlate_two_sets(cr_nchx, cr_nchx_nbc)
    # print "no Chx: correlation of barcode collapse before after {0}".format(pr)
    cr_chx, scnt_chx, dcnt_chx = get_collision_rate(cds_range, chx_dfn, chx_sfn)
    # cr_chx_nbc, scnt_chx_nbc, dcnt_chx_nbc = get_collision_rate(cds_range, chx_nbc_dfn, chx_nbc_sfn)    
    # pr = correlate_two_sets(cr_chx, cr_chx_nbc)
    # print "Chx: correlation of barcode collapse before after {0}".format(pr)
    # pr = correlate_two_sets(cr_chx_nbc, cr_nchx_nbc)
    # print "before barcode collapse: correlate chx vs nchx {0}".format(pr)

    cr_wt, scnt_wt, dcnt_wt = get_collision_rate(cds_range, wt_dfn, wt_sfn)
    cr_dom34, scnt_dom34, dcnt_dom34 = get_collision_rate(cds_range, dom34_dfn, dom34_sfn)

    # print correlation
    print "Joel: chx vs nchx (after barcode collapse)"
    print "collision rate {0:.2f}".format(correlate_two_sets(cr_chx, cr_nchx)),
    print "singlet sum: {0:.2f}".format(correlate_two_sets(scnt_chx, scnt_nchx)),
    print "doublet sum: {0:.2f}".format(correlate_two_sets(dcnt_chx, dcnt_nchx))
    print "singlet vs doublet: nchx {0:.2f}".format(correlate_two_sets(scnt_nchx, dcnt_nchx)),
    print "chx {0:.2f}".format(correlate_two_sets(scnt_chx, dcnt_chx))

    print "\nGuydosh: wt vs dom34ko"
    print "collision rate {0:.2f}".format(correlate_two_sets(cr_wt, cr_dom34)),
    print "singlet sum: {0:.2f}".format(correlate_two_sets(scnt_wt, scnt_dom34)),
    print "doublet sum: {0:.2f}".format(correlate_two_sets(dcnt_wt, dcnt_dom34))
    print "singlet vs doublet: wt {0:.2f}".format(correlate_two_sets(scnt_wt, dcnt_wt)),
    print "dom34 {0:.2f}".format(correlate_two_sets(scnt_dom34, dcnt_dom34))

    print "\nJoel chx vs Guydosh wt"
    print "collision rate {0:.2f}".format(correlate_two_sets(cr_chx, cr_wt)),
    print "singlet sum: {0:.2f}".format(correlate_two_sets(scnt_chx, scnt_wt)),
    print "doublet sum: {0:.2f}".format(correlate_two_sets(dcnt_chx, dcnt_wt))
    print "Joel nchx vs Guydosh wt"
    print "collision rate {0:.2f}".format(correlate_two_sets(cr_nchx, cr_wt)),
    print "singlet sum: {0:.2f}".format(correlate_two_sets(scnt_nchx, scnt_wt)),
    print "doublet sum: {0:.2f}".format(correlate_two_sets(dcnt_nchx, dcnt_wt))
    print "Joel chx vs Guydosh dom34"
    print "collision rate {0:.2f}".format(correlate_two_sets(cr_chx, cr_dom34)),
    print "singlet sum: {0:.2f}".format(correlate_two_sets(scnt_chx, scnt_dom34)),
    print "doublet sum: {0:.2f}".format(correlate_two_sets(dcnt_chx, dcnt_dom34))
    print "Joel nchx vs Guydosh dom34"
    print "collision rate {0:.2f}".format(correlate_two_sets(cr_nchx, cr_dom34)),
    print "singlet sum: {0:.2f}".format(correlate_two_sets(scnt_nchx, scnt_dom34)),
    print "doublet sum: {0:.2f}".format(correlate_two_sets(dcnt_nchx, dcnt_dom34))
