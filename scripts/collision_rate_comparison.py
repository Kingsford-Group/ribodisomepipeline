#!/usr/bin/env python
import os
import numpy as np
import scipy.stats
from footprint_hist_parser import parse_rlen_hist, get_cds_range
from ribomap_result_parser import parse_estimated_profile
from ribofit_utils import validate_profile, codonp_from_basep, sum_frames
from file_names import *

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['font.size'] = 12
rcParams['xtick.major.size'] = 3
rcParams['ytick.major.size'] = 3
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42

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

def get_collision_rate(cds_range, doublet_fn, singlet_fn, hc_set=False):
    print os.path.basename(singlet_fn)
    doublets = parse_rlen_hist(doublet_fn)
     # singlets (dic of ribomap outputs: rid: tid, rprofile, nprofile, mprofile)
    sp_base = parse_estimated_profile(singlet_fn)
    sp_codon = codonp_from_basep(sp_base, cds_range, sum_frames)
    # total ribo counts
    dcnt = get_total_disome_cnt(doublets, dlen_min, dlen_max, cds_range)
    print "# transcripts with doublets", len(dcnt)
    scnt = get_monocnt_from_ribomap_codon_prof(sp_codon)
    print "# transcripts with singlets", len(scnt)
    # set of transcripts to consider
    if hc_set == True:
        tid_list = high_coverage_transcript_from_ribomap_codon_prof(sp_codon)
        tcnt = len(tid_list)
        print "# transcripts with high coverage:", len(tid_list)
        tid_list = list(set(dcnt.keys())& set(tid_list))
    else:
        tid_list = list(set(dcnt.keys())& set(scnt.keys()))
        tcnt = len(scnt)
        print "# transcripts with both singlets and doublets:", len(tid_list)
    # collision rate
    cr = {}
    for tid in tid_list:
        # # no doublet, collision rate = 0
        # if tid not in dcnt:
        #     cr[tid] = 0
        #     dcnt[tid] = 0
        # # no singlet, collision rate = doublet count
        # elif tid not in scnt:
        #     cr[tid] = dcnt[tid]
        #     scnt[tid] = 0
        # else:
        cr[tid] = float(dcnt[tid])/(dcnt[tid] + scnt[tid])
    shc = { tid: scnt[tid] for tid in tid_list }
    dhc = { tid: dcnt[tid] for tid in tid_list }
    print "transcripts with non-zero doublets: {0} ({1:.0%})".format(len(cr), len(cr)/float(tcnt))
    print "average collision rate: {0:.2%}\n".format(np.mean(cr.values()))
    return cr, shc, dhc

def correlate_two_sets(s1, s2, xlabel, ylabel, title, figname):
    tid_list = set(s1.keys())&set(s2.keys())
    s1_common = [ s1[tid] for tid in tid_list ]
    s2_common = [ s2[tid] for tid in tid_list ]
    pr = scipy.stats.pearsonr(s1_common, s2_common)[0]
    sr = scipy.stats.spearmanr(s1_common, s2_common)[0]
    t = "{0}\ntotal points {1}\npearson {2:.2f} spearman {3:.2f}".format(title, len(tid_list), pr, sr)
    plt.figure()
    plt.loglog(s1_common, s2_common, c='b', marker='+', linestyle='None', alpha=0.3)
    # linear_fit = np.polyfit(s1_common, s2_common, 1)
    # fit_fn = np.poly1d(linear_fit)
    # xline =np.logspace(min(s1_common), max(s1_common), 100)
    # plt.loglog(xline, fit_fn(xline), c='r', alpha=0.7)
    plt.title(t,fontsize=12)
    plt.xlabel("{0}\ntotal points {1}".format(xlabel, len(s1)))
    plt.ylabel("{0}\ntotal points {1}".format(ylabel, len(s2)))
    plt.savefig(figname, bbox_inches='tight')
    plt.close()
    return pr

def plot_hist(nums, xlabel, title, figname):
    xmin = min(nums)
    xmax = np.percentile(nums, 95)
    plt.figure()
    ns, bins, patches = plt.hist(nums, 50, range=(xmin,xmax), color='b', alpha=0.5)
    plt.axvline(np.median(nums), c='r', lw=2, alpha=0.8)
    plt.xlabel(xlabel)
    plt.title("{0}\ntotal points {1} [{4:.2e}, {5:.2e}]\nmedian {2:.2e} std: {3:.2e}".format(title, len(nums), np.median(nums), np.std(nums), np.min(nums), np.max(nums)), fontsize=12)
    plt.savefig(figname, bbox_inches='tight')
    plt.close()

if __name__ == "__main__":
    cds_range = get_cds_range(cds_txt)
    cr_nchx, scnt_nchx, dcnt_nchx = get_collision_rate(cds_range, nchx_dfn, nchx_sfn)
    cr_chx, scnt_chx, dcnt_chx = get_collision_rate(cds_range, chx_dfn, chx_sfn)
    cr_wt, scnt_wt, dcnt_wt = get_collision_rate(cds_range, wt_dfn, wt_sfn)
    cr_dom34, scnt_dom34, dcnt_dom34 = get_collision_rate(cds_range, dom34_dfn, dom34_sfn)

    plot_hist(cr_nchx.values(), 'collision rate', 'no CHX (all)', 'cr_nchx_all_hist.pdf')
    plot_hist(scnt_nchx.values(), 'singlet loads', 'no CHX (all)', 'ss_nchx_all_hist.pdf')
    plot_hist(dcnt_nchx.values(), 'doublet loads', 'no CHX (all)', 'ds_nchx_all_hist.pdf')
    plot_hist(cr_chx.values(), 'collision rate', 'CHX (all)', 'cr_chx_all_hist.pdf')
    plot_hist(scnt_chx.values(), 'singlet loads', 'CHX (all)', 'ss_chx_all_hist.pdf')
    plot_hist(dcnt_chx.values(), 'doublet loads', 'CHX (all)', 'ds_chx_all_hist.pdf')
    plot_hist(cr_wt.values(), 'collision rate', 'wild type (all)', 'cr_wt_all_hist.pdf')
    plot_hist(scnt_wt.values(), 'singlet loads', 'wild type (all)', 'ss_wt_all_hist.pdf')
    plot_hist(dcnt_wt.values(), 'doublet loads', 'wild type (all)', 'ds_wt_all_hist.pdf')
    plot_hist(cr_dom34.values(), 'collision rate', 'DOM34 knock out (all)', 'cr_dom34_all_hist.pdf')
    plot_hist(scnt_dom34.values(), 'singlet loads', 'DOM34 knock out (all)', 'ss_dom34_all_hist.pdf')
    plot_hist(dcnt_dom34.values(), 'doublet loads', 'DOM34 knock out (all)', 'ds_dom34_all_hist.pdf')
    
    print "Joel: chx vs nchx (after barcode collapse)"
    print "collision rate {0:.2f}".format(correlate_two_sets(cr_chx, cr_nchx, 'CHX', 'no CHX', 'collision rate (all)', 'cr_joel_all.pdf')),
    print "singlet sum: {0:.2f}".format(correlate_two_sets(scnt_chx, scnt_nchx, 'CHX', 'no CHX', 'ribosome loads (all)', 'ss_joel_all.pdf')),
    print "doublet sum: {0:.2f}".format(correlate_two_sets(dcnt_chx, dcnt_nchx, 'CHX', 'no CHX', 'doublet loads (all)', 'ds_joel_all.pdf'))
    print "singlet vs doublet: nchx {0:.2f}".format(correlate_two_sets(scnt_nchx, dcnt_nchx, 'singlet loads', 'doublet loads', 'no CHX (all)', 'nchx_all.pdf')),
    print "chx {0:.2f}".format(correlate_two_sets(scnt_chx, dcnt_chx, 'singlet loads', 'doublet loads', 'CHX (all)', 'chx_all.pdf'))

    print "\nGuydosh: wt vs dom34ko"
    print "collision rate {0:.2f}".format(correlate_two_sets(cr_wt, cr_dom34, 'wild type', 'DOM34 knock out', 'collision rate (all)', 'cr_guydosh_all.pdf')),
    print "singlet sum: {0:.2f}".format(correlate_two_sets(scnt_wt, scnt_dom34, 'wild type', 'DOM34 knock out', 'ribosome loads (all)', 'ss_guydosh_all.pdf')),
    print "doublet sum: {0:.2f}".format(correlate_two_sets(dcnt_wt, dcnt_dom34, 'wild type', 'DOM34 knock out', 'doublet loads (all)', 'ds_guydosh_all.pdf'))
    print "singlet vs doublet: wt {0:.2f}".format(correlate_two_sets(scnt_wt, dcnt_wt, 'singlet loads', 'doublet loads', 'wild type (all)', 'wt_all.pdf')),
    print "dom34 {0:.2f}".format(correlate_two_sets(scnt_dom34, dcnt_dom34, 'singlet loads', 'doublet loads', 'DOM34 knock out (all)', 'dom34_all.pdf'))

    print "\nJoel chx vs Guydosh wt"
    print "collision rate {0:.2f}".format(correlate_two_sets(cr_chx, cr_wt, 'CHX', 'wild type', 'collision rate (all)', 'cr_chx_wt_all.pdf')),
    print "singlet sum: {0:.2f}".format(correlate_two_sets(scnt_chx, scnt_wt, 'CHX', 'wild type', 'ribosome loads (all)', 'ss_chx_wt_all.pdf')),
    print "doublet sum: {0:.2f}".format(correlate_two_sets(dcnt_chx, dcnt_wt, 'CHX', 'wild type', 'doublet loads (all)', 'ds_chx_wt_all.pdf'))
    print "Joel nchx vs Guydosh wt"
    print "collision rate {0:.2f}".format(correlate_two_sets(cr_nchx, cr_wt, 'no CHX', 'wild type', 'collision rate (all)', 'cr_nchx_wt_all.pdf')),
    print "singlet sum: {0:.2f}".format(correlate_two_sets(scnt_nchx, scnt_wt, 'no CHX', 'wild type', 'ribosome loads (all)', 'ss_nchx_wt_all.pdf')),
    print "doublet sum: {0:.2f}".format(correlate_two_sets(dcnt_nchx, dcnt_wt, 'no CHX', 'wild type', 'doublet loads (all)', 'ds_nchx_wt_all.pdf'))
    print "Joel chx vs Guydosh dom34"
    print "collision rate {0:.2f}".format(correlate_two_sets(cr_chx, cr_dom34, 'CHX', 'DOM34 knock out', 'collision rate (all)', 'cr_chx_dom34_all.pdf')),
    print "singlet sum: {0:.2f}".format(correlate_two_sets(scnt_chx, scnt_dom34, 'CHX', 'DOM34 knock out', 'ribosome loads (all)', 'ss_chx_dom34_all.pdf')),
    print "doublet sum: {0:.2f}".format(correlate_two_sets(dcnt_chx, dcnt_dom34, 'CHX', 'DOM34 knock out', 'doublet loads (all)', 'ds_chx_dom34_all.pdf'))
    print "Joel nchx vs Guydosh dom34"
    print "collision rate {0:.2f}".format(correlate_two_sets(cr_nchx, cr_dom34, 'no CHX', 'DOM34 knock out', 'collision rate (all)', 'cr_nchx_dom34_all.pdf')),
    print "singlet sum: {0:.2f}".format(correlate_two_sets(scnt_nchx, scnt_dom34, 'no CHX', 'DOM34 knock out', 'ribosome loads (all)', 'ss_nchx_dom34_all.pdf')),
    print "doublet sum: {0:.2f}".format(correlate_two_sets(dcnt_nchx, dcnt_dom34, 'no CHX', 'DOM34 knock out', 'doublet loads (all)', 'ds_nchx_dom34_all.pdf'))

    cr_nchx, scnt_nchx, dcnt_nchx = get_collision_rate(cds_range, nchx_dfn, nchx_sfn, True)
    cr_chx, scnt_chx, dcnt_chx = get_collision_rate(cds_range, chx_dfn, chx_sfn, True)
    cr_wt, scnt_wt, dcnt_wt = get_collision_rate(cds_range, wt_dfn, wt_sfn, True)
    cr_dom34, scnt_dom34, dcnt_dom34 = get_collision_rate(cds_range, dom34_dfn, dom34_sfn, True)

    plot_hist(cr_nchx.values(), 'collision rate', 'no CHX (high coverage)', 'cr_nchx_hc_hist.pdf')
    plot_hist(scnt_nchx.values(), 'singlet loads', 'no CHX (high coverage)', 'ss_nchx_hc_hist.pdf')
    plot_hist(dcnt_nchx.values(), 'doublet loads', 'no CHX (high coverage)', 'ds_nchx_hc_hist.pdf')
    plot_hist(cr_chx.values(), 'collision rate', 'CHX (high coverage)', 'cr_chx_hc_hist.pdf')
    plot_hist(scnt_chx.values(), 'singlet loads', 'CHX (high coverage)', 'ss_chx_hc_hist.pdf')
    plot_hist(dcnt_chx.values(), 'doublet loads', 'CHX (high coverage)', 'ds_chx_hc_hist.pdf')
    plot_hist(cr_wt.values(), 'collision rate', 'wild type (high coverage)', 'cr_wt_hc_hist.pdf')
    plot_hist(scnt_wt.values(), 'singlet loads', 'wild type (high coverage)', 'ss_wt_hc_hist.pdf')
    plot_hist(dcnt_wt.values(), 'doublet loads', 'wild type (high coverage)', 'ds_wt_hc_hist.pdf')
    plot_hist(cr_dom34.values(), 'collision rate', 'DOM34 knock out (high coverage)', 'cr_dom34_hc_hist.pdf')
    plot_hist(scnt_dom34.values(), 'singlet loads', 'DOM34 knock out (high coverage)', 'ss_dom34_hc_hist.pdf')
    plot_hist(dcnt_dom34.values(), 'doublet loads', 'DOM34 knock out (high coverage)', 'ds_dom34_hc_hist.pdf')

    print "Joel: chx vs nchx (after barcode collapse)"
    print "collision rate {0:.2f}".format(correlate_two_sets(cr_chx, cr_nchx, 'CHX', 'no CHX', 'collision rate (high coverage)', 'cr_joel_hc.pdf')),
    print "singlet sum: {0:.2f}".format(correlate_two_sets(scnt_chx, scnt_nchx, 'CHX', 'no CHX', 'ribosome loads (high coverage)', 'ss_joel_hc.pdf')),
    print "doublet sum: {0:.2f}".format(correlate_two_sets(dcnt_chx, dcnt_nchx, 'CHX', 'no CHX', 'doublet loads (high coverage)', 'ds_joel_hc.pdf'))
    print "singlet vs doublet: nchx {0:.2f}".format(correlate_two_sets(scnt_nchx, dcnt_nchx, 'singlet loads', 'doublet loads', 'no CHX (high coverage)', 'nchx_hc.pdf')),
    print "chx {0:.2f}".format(correlate_two_sets(scnt_chx, dcnt_chx, 'singlet loads', 'doublet loads', 'CHX (high coverage)', 'chx_hc.pdf'))

    print "\nGuydosh: wt vs dom34ko"
    print "collision rate {0:.2f}".format(correlate_two_sets(cr_wt, cr_dom34, 'wild type', 'DOM34 knock out', 'collision rate (high coverage)', 'cr_guydosh_hc.pdf')),
    print "singlet sum: {0:.2f}".format(correlate_two_sets(scnt_wt, scnt_dom34, 'wild type', 'DOM34 knock out', 'ribosome loads (high coverage)', 'ss_guydosh_hc.pdf')),
    print "doublet sum: {0:.2f}".format(correlate_two_sets(dcnt_wt, dcnt_dom34, 'wild type', 'DOM34 knock out', 'doublet loads (high coverage)', 'ds_guydosh_hc.pdf'))
    print "singlet vs doublet: wt {0:.2f}".format(correlate_two_sets(scnt_wt, dcnt_wt, 'singlet loads', 'doublet loads', 'wild type (high coverage)', 'wt_hc.pdf')),
    print "dom34 {0:.2f}".format(correlate_two_sets(scnt_dom34, dcnt_dom34, 'singlet loads', 'doublet loads', 'DOM34 knock out (high coverage)', 'dom34_hc.pdf'))

    print "\nJoel chx vs Guydosh wt"
    print "collision rate {0:.2f}".format(correlate_two_sets(cr_chx, cr_wt, 'CHX', 'wild type', 'collision rate (high coverage)', 'cr_chx_wt_hc.pdf')),
    print "singlet sum: {0:.2f}".format(correlate_two_sets(scnt_chx, scnt_wt, 'CHX', 'wild type', 'ribosome loads (high coverage)', 'ss_chx_wt_hc.pdf')),
    print "doublet sum: {0:.2f}".format(correlate_two_sets(dcnt_chx, dcnt_wt, 'CHX', 'wild type', 'doublet loads (high coverage)', 'ds_chx_wt_hc.pdf'))
    print "Joel nchx vs Guydosh wt"
    print "collision rate {0:.2f}".format(correlate_two_sets(cr_nchx, cr_wt, 'no CHX', 'wild type', 'collision rate (high coverage)', 'cr_nchx_wt_hc.pdf')),
    print "singlet sum: {0:.2f}".format(correlate_two_sets(scnt_nchx, scnt_wt, 'no CHX', 'wild type', 'ribosome loads (high coverage)', 'ss_nchx_wt_hc.pdf')),
    print "doublet sum: {0:.2f}".format(correlate_two_sets(dcnt_nchx, dcnt_wt, 'no CHX', 'wild type', 'doublet loads (high coverage)', 'ds_nchx_wt_hc.pdf'))
    print "Joel chx vs Guydosh dom34"
    print "collision rate {0:.2f}".format(correlate_two_sets(cr_chx, cr_dom34, 'CHX', 'DOM34 knock out', 'collision rate (high coverage)', 'cr_chx_dom34_hc.pdf')),
    print "singlet sum: {0:.2f}".format(correlate_two_sets(scnt_chx, scnt_dom34, 'CHX', 'DOM34 knock out', 'ribosome loads (high coverage)', 'ss_chx_dom34_hc.pdf')),
    print "doublet sum: {0:.2f}".format(correlate_two_sets(dcnt_chx, dcnt_dom34, 'CHX', 'DOM34 knock out', 'doublet loads (high coverage)', 'ds_chx_dom34_hc.pdf'))
    print "Joel nchx vs Guydosh dom34"
    print "collision rate {0:.2f}".format(correlate_two_sets(cr_nchx, cr_dom34, 'no CHX', 'DOM34 knock out', 'collision rate (high coverage)', 'cr_nchx_dom34_hc.pdf')),
    print "singlet sum: {0:.2f}".format(correlate_two_sets(scnt_nchx, scnt_dom34, 'no CHX', 'DOM34 knock out', 'ribosome loads (high coveage)', 'ss_nchx_dom34_hc.pdf')),
    print "doublet sum: {0:.2f}".format(correlate_two_sets(dcnt_nchx, dcnt_dom34, 'no CHX', 'DOM34 knock out', 'doublet loads (high coverage)', 'ds_nchx_dom34_hc.pdf'))



    # cr_nchx_nbc, scnt_nchx_nbc, dcnt_nchx_nbc = get_collision_rate(cds_range, nchx_nbc_dfn, nchx_nbc_sfn)    
    # pr = correlate_two_sets(cr_nchx, cr_nchx_nbc)
    # print "no Chx: correlation of barcode collapse before after {0}".format(pr)
    # cr_chx_nbc, scnt_chx_nbc, dcnt_chx_nbc = get_collision_rate(cds_range, chx_nbc_dfn, chx_nbc_sfn)    
    # pr = correlate_two_sets(cr_chx, cr_chx_nbc)
    # print "Chx: correlation of barcode collapse before after {0}".format(pr)
    # pr = correlate_two_sets(cr_chx_nbc, cr_nchx_nbc)
    # print "before barcode collapse: correlate chx vs nchx {0}".format(pr)
