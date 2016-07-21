#!/usr/bin/env python
import numpy as np
import scipy.stats
from doublet_profile import generate_codon_profile_from_rlen_hist, get_tid2codonp_from_ribomap_base, get_codon_profile_from_deblur_profile
from significant_doublet_count import get_window_cnt, read_pvals_from_file
from io_utils import get_cds_range, get_tseq
from file_names import *

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['font.size'] = 20
rcParams['xtick.major.size'] = 5
rcParams['ytick.major.size'] = 5

#==============================
# utils
#==============================
def get_significant_peaks(tid2pvals, sig_cutoff=0.1):
    tot_cnt = 0
    sig_cnt = 0
    tid2sig = {}
    tid2rest = {}
    for tid in tid2pvals:
        for pos, pval in tid2pvals[tid]:
            tot_cnt += 1
            if pval <= sig_cutoff:
                tid2sig.setdefault(tid, []).append(pos)
                sig_cnt += 1
            else:
                tid2rest.setdefault(tid, []).append(pos)
    print "total peaks: {0} significant peaks {1} ({2:.2%})".format(tot_cnt, sig_cnt, float(sig_cnt)/tot_cnt)
    return tid2sig, tid2rest

def get_sum(pos, prof):
    return np.sum(prof)

def get_mean(pos, prof):
    return np.mean(prof)

def get_mean_adjusted(pos, prof):
    return prof[pos]/np.mean(prof)

def get_relative_pos(pos, prof):
    return float(pos)/len(prof)

def max_upstream(pos, prof):
    iend = pos - sd_distance
    if iend < 0 : return None
    istart = iend - window_size
    if istart >= len(prof): return None
    if istart < 0: istart = 0
    if istart == iend: return 0
    return np.max(prof[istart:iend])

def max_downstream(pos, prof):
    iend = pos - ds_distance
    if iend < 0 : return None
    istart = iend - window_size
    if istart >= len(prof): return None
    if istart < 0: istart = 0
    if istart == iend: return 0
    return np.max(prof[istart:iend])

def max_upstream_suprise(pos, prof):
    iend = pos - sd_distance
    if iend < 0 : return None
    istart = iend - window_size
    if istart >= len(prof): return None
    if istart < 0: istart = 0
    if istart == iend: return 0
    return np.max(prof[istart:iend])/np.mean(prof)

def max_downstream_suprise(pos, prof):
    iend = pos - ds_distance
    if iend < 0 : return None
    istart = iend - window_size
    if istart >= len(prof): return None
    if istart < 0: istart = 0
    if istart == iend: return 0
    return np.max(prof[istart:iend])/np.mean(prof)

def get_values(tid2pos, profile_list, func=None):
    vec = []
    for tid in tid2pos:
        for pos in tid2pos[tid]:
            if func == None:
                vec.append(profile_list[tid][pos])
            else:
                vec.append(func(pos, profile_list[tid]))
    return vec

def get_collision_rates(tid2pos, sprof, dprof, peak_type):
    vec = []
    for tid in tid2pos:
        for pos in tid2pos[tid]:
            if peak_type == 'singlet':
                dcnt = max_upstream(pos, dprof[tid])
                scnt = sprof[tid][pos]
            elif peak_type == 'doublet':
                dcnt = dprof[tid][pos]
                scnt = max_downstream(pos, sprof[tid])
            else:
                print "peak type {0} not supported!"
                exit(1)
            if scnt != 0:
                cr = float(dcnt)/scnt
                vec.append(cr)
    return vec

#==============================
# jam vs non-jam
#==============================
def boxplot_compare(vec_sig, vec_rest, ylabel, figname):
    plt.figure()
    plt.boxplot([ vec_sig, vec_rest ], whis='range', labels=['peaks with \nsignificant collision', 'rest' ])
    up_val = max(np.percentile(vec_sig, 80), np.percentile(vec_rest, 80))
    low_val = np.floor(min(np.percentile(vec_sig, 0), np.percentile(vec_rest, 0)))
    plt.ylabel(ylabel)
    plt.ylim((low_val,up_val))
    plt.savefig(figname, bbox_inches='tight')
    plt.close()

def significance_analysis_pipeline(pfname, dfname, sfname, oprefix):
    """ two group comparision: jam vs nonjam"""
    tid2sig, tid2rest = get_sig_peaks_from_file(pfname)
    print "getting doublet profile..."
    dcp = generate_codon_profile_from_rlen_hist(dfname, cds_range)
    print "getting singlet profile..."
    scp = get_tid2codonp_from_ribomap_base(sfname, cds_range)
    vec_sig = get_values(tid2sig, scp)
    vec_rest = get_values(tid2rest, scp)
    boxplot_compare(vec_sig, vec_rest, 'singlet peak', oprefix+"_sp.pdf")
    vec_sig = get_values(tid2sig, scp, get_mean)
    vec_rest = get_values(tid2rest, scp, get_mean)
    boxplot_compare(vec_sig, vec_rest, 'singlet coverage', oprefix+"_sm.pdf")
    vec_sig = get_values(tid2sig, scp, get_sum)
    vec_rest = get_values(tid2rest, scp, get_sum)
    boxplot_compare(vec_sig, vec_rest, 'singlet load', oprefix+"_ss.pdf")
    vec_sig = get_values(tid2sig, scp, get_mean_adjusted)
    vec_rest = get_values(tid2rest, scp, get_mean_adjusted)
    boxplot_compare(vec_sig, vec_rest, 'singlet suprise', oprefix+"_spadj.pdf")
    vec_sig = get_values(tid2sig, scp, get_mean)
    vec_rest = get_values(tid2rest, scp, get_mean)
    boxplot_compare(vec_sig, vec_rest, 'singlet coverage', oprefix+"_sm.pdf")
    vec_sig = get_values(tid2sig, scp, get_relative_pos)
    vec_rest = get_values(tid2rest, scp, get_relative_pos)
    boxplot_compare(vec_sig, vec_rest, 'relative peak position', oprefix+"_srpos.pdf")
    vec_sig = get_values(tid2sig, dcp, get_window_cnt)
    vec_rest = get_values(tid2rest, dcp, get_window_cnt)
    boxplot_compare(vec_sig, vec_rest, 'doublet peak', oprefix+"_dp.pdf")
    vec_sig = get_values(tid2sig, dcp, get_mean)
    vec_rest = get_values(tid2rest, dcp, get_mean)
    boxplot_compare(vec_sig, vec_rest, 'doublet coverage', oprefix+"_dm.pdf")
    vec_sig = get_values(tid2sig, dcp, get_sum)
    vec_rest = get_values(tid2rest, dcp, get_sum)
    boxplot_compare(vec_sig, vec_rest, 'doublet load', oprefix+"_ds.pdf")
    cr_sig = get_collision_rates(tid2sig, scp, dcp)
    cr_rest = get_collision_rates(tid2rest, scp, dcp)
    boxplot_compare(cr_sig, cr_rest, 'collision rate', oprefix+"_cr.pdf")

#==============================
# code sanity check
#==============================
def significance_cmp_before_after(p1fname, p2fname):
    tid2pvals1 = read_pvals_from_file(p1fname)
    tid2sig1, tid2rest = get_significant_peaks(tid2pvals1)
    tid2pvals2 = read_pvals_from_file(p2fname)
    tid2sig2, tid2rest = get_significant_peaks(tid2pvals2)
    for tid in tid2sig1:
        if tid not in tid2sig2:
            print "oops tid {0} not exists in file {1}!".format(tid, p2fname)
            print tid2sig1[tid]
            print tid2pvals1[tid]
            print tid2pvals2[tid]
            continue
        if tid2sig1[tid] != tid2sig2[tid]:
            print "oops two files disagree!"
            print tid, tid2sig1[tid], tid2sig2[tid]
            print tid2pvals1[tid]
            print tid2pvals2[tid]

#================================
# compare joint jam, jam, non-jam
#================================
def overlapped_peaks(dlist, slist, distance, window_size, peak_type):
    if peak_type == 'singlet':
        convert_list = dlist
        ref_list = slist
    elif peak_type == 'doublet':
        convert_list = slist
        ref_list = dlist
    else:
        print "unsupported peak type {0}!".format(peak_type)
        exit(1)
    cmp_list = set([])
    for pos in convert_list:
        for i in xrange(pos-distance-window_size, pos-distance):
            cmp_list.add(i)
    return list(set(ref_list) & cmp_list)
            
def batch_overlapped_peaks(tid2speak, tid2dpeak, distance, window_size, peak_type):
    double_sig = {}
    scnt = sum(map(len, tid2speak.values()))
    dcnt = sum(map(len, tid2dpeak.values()))
    for tid in tid2dpeak: 
        if tid not in tid2speak: continue
        op = overlapped_peaks(tid2dpeak[tid], tid2speak[tid], distance, window_size, peak_type)
        if len(op)!=0: double_sig[tid] = op
    ccnt = sum(map(len, double_sig.values()))
    print "significant singlet {0} doublet {1} joint {2} ({3:.0%})".format(scnt, dcnt, ccnt, float(ccnt)/dcnt)
    return double_sig

def get_sig_peaks_from_file(pfname, cutoff=0.1):
    tid2pvals = read_pvals_from_file(pfname)
    tid2sig, tid2rest = get_significant_peaks(tid2pvals, cutoff)
    return tid2sig, tid2rest

def prepare_peak_list(psfname, pdfname, distance, window_size, peak_type):
    """
    separate peaks into 5 groups:
    singlet-jam, singlet-nonjam, doublet-jam, doublet-non-jam, joint-jam
    """
    print "getting significant singlet peaks..."
    tid2ssig, tid2srest = get_sig_peaks_from_file(psfname)
    print "getting significant doublet peaks..."
    tid2dsig, tid2drest = get_sig_peaks_from_file(pdfname)
    print "getting joint significant peaks..."
    tid2jsig = batch_overlapped_peaks(tid2ssig, tid2dsig, distance, window_size, peak_type)
    return tid2ssig, tid2srest, tid2dsig, tid2drest, tid2jsig

def group_boxplots(group_data, group_label, ylabel, figname):
    plt.figure()
    plt.boxplot(group_data, whis='range', labels=group_label)
    up_val = max(map(lambda v: np.percentile(v, 80), group_data))
    low_val = np.floor(min(map(lambda v: np.percentile(v, 0), group_data)))
    plt.ylabel(ylabel)
    plt.ylim((low_val,up_val))
    plt.savefig(figname, bbox_inches='tight')
    plt.close()

def plot_dspeaks(dprof, sprof, dpeak, speak, dpsig, spsig, tid, figname, grid=False):
    cds_len = len(dprof)
    assert len(dprof) == len(sprof)
    x_pos = np.arange(cds_len)
    xgrid = np.arange(0,cds_len, 3)
    width = min(300, int(cds_len/10.0*1.5))
    plt.figure(figsize=(width, 15))
    plt.subplot(2,1,1)
    if grid == True:
        plt.vlines(xgrid, [0]*len(xgrid), [max(dprof)+1]*len(xgrid), 'r', linewidths=5, alpha=0.1)
    plt.vlines(x_pos, [0]*cds_len, dprof, 'k', linewidths=5, alpha=0.2)
    plt.xticks(np.arange(0,cds_len, 10))
    ypeaks = [ dprof[i] for i in dpeak ]
    plt.scatter(dpeak, ypeaks, s=20, c='b', marker="D", edgecolors='face', alpha=0.5)
    ypeaks = [ dprof[i] for i in dpsig ]
    plt.scatter(dpsig, ypeaks, s=50, c='r', marker="D", edgecolors='face')
    plt.xlim((-1, cds_len))
    plt.ylim((0,max(dprof)+1))
    plt.title("{0}\ndoublet profile".format(tid))
    plt.subplot(2,1,2)
    if grid == True:
        plt.vlines(xgrid, [0]*len(xgrid), [max(sprof)+1]*len(xgrid), 'r', linewidths=5, alpha=0.1)
    plt.vlines(x_pos, [0]*cds_len, sprof, 'k', linewidths=5, alpha=0.2)
    plt.xticks(np.arange(0,cds_len, 10))
    ypeaks = [ sprof[i] for i in speak ]
    plt.scatter(speak, ypeaks, s=20, c='b', marker="D", edgecolors='face', alpha=0.5)
    ypeaks = [ sprof[i] for i in spsig ]
    plt.scatter(spsig, ypeaks, s=50, c='r', marker="D", edgecolors='face')
    plt.xlim((-1, cds_len))
    plt.ylim((0, max(sprof)+1))
    cr = float(sum(dprof))/(sum(dprof)+sum(sprof))
    plt.title("collision rate {0:.2%}\nsinglet profile".format(cr))
    plt.xlabel("position")
    try:
        plt.savefig(figname, bbox_inches='tight')
    except ValueError:
        print "render error occurs! try without tight layout"
        try:
            plt.savefig(figname, bbox_inches='tight')
        except ValueError:
            print "render error occurs without tight layout!"
    plt.close()
    print "{0} len: {1} average doublet: {2:.2f} collision rate: {3:.2%}".format(tid,len(sprof), np.mean(dprof), cr)

def double_significance_pipeline(pdfname, psfname, distance, window_size, peak_type, dfname, sfname, oprefix, multimap):
    print "significant singlet"
    tid2ssig, tid2srest = get_sig_peaks_from_file(psfname)
    print "significant doublet"
    tid2dsig, tid2drest = get_sig_peaks_from_file(pdfname)
    tid2jsig = batch_overlapped_peaks(tid2ssig, tid2dsig, distance, window_size, peak_type)
    print "getting doublet profile..."
    dcp = generate_codon_profile_from_rlen_hist(dfname, cds_range)
    print "getting singlet profile..."
    if multimap == True:
        scp = get_tid2codonp_from_ribomap_base(sfname, cds_range)
    else:
        scp = get_codon_profile_from_deblur_profile(sfname)
    print "boxplot comparision..."
    group_label = ['significant\njoint\npeaks','significant\nsinglet\npeaks', 'insignificant\nsinglet\npeaks', 'significant\ndoublet\npeaks', 'insignificant\ndoublet\npeaks']
    tid2sosig = {}
    for tid in tid2ssig:
        if tid not in tid2jsig: tid2sosig[tid] = tid2ssig[tid]
        else: tid2sosig[tid] = list(set(tid2ssig[tid])-set(tid2jsig[tid]))
    jsig = get_values(tid2jsig, scp)
    ssig = get_values(tid2sosig, scp)
    srest = get_values(tid2srest, scp)
    dsig = get_values(tid2dsig, dcp)
    drest = get_values(tid2drest, dcp)
    group_data = [ jsig, ssig, srest, dsig, drest ]
    group_boxplots(group_data, group_label, 'peak values', oprefix+"_peak.pdf")
    jsig = get_values(tid2jsig, scp, get_mean)
    ssig = get_values(tid2sosig, scp, get_mean)
    srest = get_values(tid2srest, scp, get_mean)
    dsig = get_values(tid2dsig, dcp, get_mean)
    drest = get_values(tid2drest, dcp, get_mean)
    group_data = [ jsig, ssig, srest, dsig, drest ]
    group_boxplots(group_data, group_label, 'profile coverage', oprefix+"_cov.pdf")
    jsig = get_values(tid2jsig, scp, get_sum)
    ssig = get_values(tid2sosig, scp, get_sum)
    srest = get_values(tid2srest, scp, get_sum)
    dsig = get_values(tid2dsig, dcp, get_sum)
    drest = get_values(tid2drest, dcp, get_sum)
    group_data = [ jsig, ssig, srest, dsig, drest ]
    group_boxplots(group_data, group_label, 'footprint loads', oprefix+"_lds.pdf")
    jsig = get_values(tid2jsig, scp, get_mean_adjusted)
    ssig = get_values(tid2sosig, scp, get_mean_adjusted)
    srest = get_values(tid2srest, scp, get_mean_adjusted)
    dsig = get_values(tid2dsig, dcp, get_mean_adjusted)
    drest = get_values(tid2drest, dcp, get_mean_adjusted)
    group_data = [ jsig, ssig, srest, dsig, drest ]
    group_boxplots(group_data, group_label, 'peak suprise', oprefix+"_padj.pdf")
    jsig = get_values(tid2jsig, scp, get_relative_pos)
    ssig = get_values(tid2sosig, scp, get_relative_pos)
    srest = get_values(tid2srest, scp, get_relative_pos)
    dsig = get_values(tid2dsig, dcp, get_relative_pos)
    drest = get_values(tid2drest, dcp, get_relative_pos)
    group_data = [ jsig, ssig, srest, dsig, drest ]
    group_boxplots(group_data, group_label, 'relative peak position', oprefix+"_rpos.pdf")

    jsig = get_values(tid2jsig, dcp, max_upstream_suprise)
    ssig = get_values(tid2sosig, dcp, max_upstream_suprise)
    srest = get_values(tid2srest, dcp, max_upstream_suprise)
    dsig = get_values(tid2dsig, scp, max_downstream_suprise)
    drest = get_values(tid2drest, scp, max_downstream_suprise)
    group_data = [ jsig, ssig, srest, dsig, drest ]
    group_boxplots(group_data, group_label, 'max up/downstream suprise', oprefix+"_maxpadj.pdf")

    jsig = get_collision_rates(tid2jsig, scp, dcp, 'singlet')
    ssig = get_collision_rates(tid2sosig, scp, dcp, 'singlet')
    srest = get_collision_rates(tid2srest, scp, dcp, 'singlet')
    dsig = get_collision_rates(tid2dsig, scp, dcp, 'doublet')
    drest = get_collision_rates(tid2drest, scp, dcp, 'doublet')
    group_data = [ jsig, ssig, srest, dsig, drest ]
    group_boxplots(group_data, group_label, 'collision rate', oprefix+"_cr.pdf")

    # print "joint tids: ", len(tid2jsig)
    # for tid in tid2jsig:
    #     print "doublet pvals:", tid2dpvals[tid]
    #     print "singlet pvals:", tid2spvals[tid]
    #     print "joint significant peaks:", tid2jsig[tid]
    #     dpeak = [ pos for pos,_ in tid2dpvals[tid] ]
    #     speak = [ pos for pos,_ in tid2spvals[tid] ]
    #     plot_dspeaks(dcp[tid], scp[tid], dpeak, speak, tid2dsig[tid], tid2ssig[tid], tid, "{0}.pdf".format(tid))
    #     break

#==============================
# jam reproducibility
#==============================
def closest_peak(p1, p2):
    return [ min(np.abs(pos-np.array(p2))) for pos in p1 ]
    
def compare_peak_list(tid2p1, tid2p2, fn1, fn2, figname, close_range=10):
    tid_list = list(set(tid2p1.keys())&set(tid2p2.keys()))
    print "overlapping transcripts: {0}".format(len(tid_list))
    plt.figure(figsize=(12,6))
    plt.subplot(121)
    dlist = []
    for tid in tid_list:
        dlist.extend(closest_peak(tid2p1[tid], tid2p2[tid]))
    close_cnt = np.sum(np.array(dlist)<=close_range)
    ns, bins, patches = plt.hist(dlist, bins=range(0,101,5), color='b', alpha=0.5)
    plt.xlabel("distance to closest peak")
    plt.ylabel("# peaks")
    plt.title("{0} compare to {1}\ntotal peaks in {0}: {2} [{3}, {4}]\nmatched peak within 10 codons: {5} ({6:.0%})".format(fn1, fn2, len(dlist), min(dlist), max(dlist), close_cnt, float(close_cnt)/len(dlist)))
    plt.subplot(122)
    dlist = []
    for tid in tid_list:
        dlist.extend(closest_peak(tid2p2[tid], tid2p1[tid]))
    close_cnt = np.sum(np.array(dlist)<=close_range)
    ns, bins, patches = plt.hist(dlist, bins=range(0,101,5), color='b', alpha=0.5)
    plt.xlabel("distance to closest peak")
    plt.ylabel("# peaks")
    plt.title("{0} compare to {1}\ntotal peaks in {0}: {2} [{3}, {4}]\nmatched peak within 10 codons: {5} ({6:.0%})".format(fn2, fn1, len(dlist), min(dlist), max(dlist), close_cnt, float(close_cnt)/len(dlist)))
    plt.savefig(figname, bbox_inches='tight')
    plt.close()

def pair_peak_compare_pipeline(p1sfname, p1dfname, p2sfname, p2dfname, distance, window_size, peak_type, oprfx, fn1, fn2):
    tid2ssig1, tid2srest1, tid2dsig1, tid2drest1, tid2jsig1 = prepare_peak_list(p1sfname, p1dfname, distance, window_size, peak_type)
    tid2ssig2, tid2srest2, tid2dsig2, tid2drest2, tid2jsig2 = prepare_peak_list(p2sfname, p2dfname, distance, window_size, peak_type)
    figname = "{0}singlet_sig_{1}_{2}.pdf".format(oprfx, fn1, fn2)
    compare_peak_list(tid2ssig1, tid2ssig2, fn1, fn2, figname)
    figname = "{0}doublet_sig_{1}_{2}.pdf".format(oprfx, fn1, fn2)
    compare_peak_list(tid2dsig1, tid2dsig2, fn1, fn2, figname)
    figname = "{0}joint_sig_{1}_{2}.pdf".format(oprfx, fn1, fn2)
    compare_peak_list(tid2jsig1, tid2jsig2, fn1, fn2, figname)
    figname = "{0}singlet_insig_{1}_{2}.pdf".format(oprfx, fn1, fn2)
    compare_peak_list(tid2srest1, tid2srest2, fn1, fn2, figname)
    figname = "{0}doublet_insig_{1}_{2}.pdf".format(oprfx, fn1, fn2)
    compare_peak_list(tid2drest1, tid2drest2, fn1, fn2, figname)
    
def reproducibility_sics_nonsics(p1fname, p2fname, figname, close_range=10):
    tid2sig1, tid2rest1 = get_sig_peaks_from_file(p1fname)
    tid2sig2, tid2rest2 = get_sig_peaks_from_file(p2fname)
    
    "print significant set"
    tid_list = list(set(tid2sig1.keys())&set(tid2sig2.keys()))
    print "overlapping transcripts: {0}".format(len(tid_list))
    plt.figure(figsize=(12,6))
    plt.subplot(121)
    dlist_sig = []
    for tid in tid_list:
        dlist_sig.extend(closest_peak(tid2sig1[tid], tid2sig2[tid]))
    close_cnt = np.sum(np.array(dlist_sig)<=close_range)
    ns, bins, patches = plt.hist(dlist_sig, bins=range(0,101,5), color='b', alpha=0.5)
    plt.xlabel("distance to closest peak")
    plt.ylabel("# peaks")
    plt.title("SICS",size=20)
    print "total peaks file 1: {0} [{1}, {2}]".format(len(dlist_sig), min(dlist_sig), max(dlist_sig))
    print "matched peak within 10 codons: {0} ({1:.0%})".format(close_cnt, float(close_cnt)/len(dlist_sig))
    "print in-significant set"
    tid_list = list(set(tid2rest1.keys())&set(tid2rest2.keys()))
    plt.subplot(122)
    dlist_insig = []
    for tid in tid_list:
        dlist_insig.extend(closest_peak(tid2rest1[tid], tid2rest2[tid]))
    close_cnt = np.sum(np.array(dlist_insig)<=close_range)
    ns, bins, patches = plt.hist(dlist_insig, bins=range(0,101,5), color='b', alpha=0.5)
    plt.xlabel("distance to closest peak")
    plt.ylabel("# peaks")
    plt.title("nonSICS",size=20)
    print "total peaks in file 1: {0} [{1}, {2}]".format(len(dlist_insig), min(dlist_insig), max(dlist_insig))
    print "matched peak within 10 codons: {0} ({1:.0%})".format(close_cnt, float(close_cnt)/len(dlist_insig))
    plt.tight_layout()
    plt.savefig(figname, bbox_inches='tight')
    plt.close()
    print "SICS VS nonSICS: mwu p-value:", scipy.stats.mannwhitneyu(dlist_sig, dlist_insig)

if __name__ == "__main__":
    cds_range = get_cds_range(cds_txt) 

    figname = figure_dir+"sics_match_nchx_wt.pdf"
    reproducibility_sics_nonsics(nchx_psfn, wt_psfn, figname)

    """
    distance = sd_distance
    peak_type = 'doublet'

    # reproducibility on jam locations
    pair_peak_compare_pipeline(nchx_psfn, nchx_pdfn, chx_psfn, chx_pdfn, distance, window_size, peak_type, figure_dir, 'nchx', 'chx')
    pair_peak_compare_pipeline(nchx_psfn, nchx_pdfn, wt_psfn, wt_pdfn, distance, window_size, peak_type, figure_dir, 'nchx', 'wt')
    pair_peak_compare_pipeline(nchx_psfn, nchx_pdfn, dom34_psfn, dom34_pdfn, distance, window_size, peak_type, figure_dir, 'nchx', 'dom34')
    pair_peak_compare_pipeline(chx_psfn, chx_pdfn, wt_psfn, wt_pdfn, distance, window_size, peak_type, figure_dir, 'chx', 'wt')
    pair_peak_compare_pipeline(chx_psfn, chx_pdfn, dom34_psfn, dom34_pdfn, distance, window_size, peak_type, figure_dir, 'chx', 'dom34')
    pair_peak_compare_pipeline(wt_psfn, wt_pdfn, dom34_psfn, dom34_pdfn, distance, window_size, peak_type, figure_dir, 'wt', 'dom34')

    # jam features
    double_significance_pipeline(nchx_pdfn, nchx_psfn, distance, window_size, peak_type, nchx_dfn, nchx_sfn, figure_dir+"nchx_joint", multimap)
    double_significance_pipeline(chx_pdfn, chx_psfn, distance, window_size, peak_type, chx_dfn, chx_sfn, figure_dir+"chx_joint", multimap)
    double_significance_pipeline(wt_pdfn, wt_psfn, distance, window_size, peak_type, wt_dfn, wt_sfn, figure_dir+"wt_joint", multimap)
    double_significance_pipeline(dom34_pdfn, dom34_psfn, distance, window_size, peak_type, dom34_dfn, dom34_sfn, figure_dir+"dom34_joint", multimap)

    # eyeballing sanity check to see whether pipline is broken
    significance_cmp_before_after("../ds_cmp/nchx_pval.txt", "../ds_cmp/nchx_singlet_pval.txt")
    significance_analysis_pipeline(nchx_dfn, nchx_sfn, "../ds_cmp/nchx_pval")
    significance_analysis_pipeline(chx_pfn, chx_dfn, chx_sfn, "../ds_cmp/chx_pval")
    significance_analysis_pipeline(wt_pfn, wt_dfn, wt_sfn, "../ds_cmp/wt_pval")
    significance_analysis_pipeline(dom34_pfn, dom34_dfn, dom34_sfn, "../ds_cmp/dom34_pval")
    """
