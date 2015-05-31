#!/usr/bin/env python
import sys
import numpy as np
import scipy.stats
from codon_table import *
from ribomap_result_parser import *
from peak_cluster import threshold, identify_peaks, validate_profile
from ribofit_utils import codonp_from_basep, sum_frames

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['font.size'] = 12
rcParams['xtick.major.size'] = 5
rcParams['ytick.major.size'] = 5
cmap = matplotlib.cm.Paired

def get_peaks_from_profile(prof, start, stop, valid_func, peak_func):
    print "calling peaks..."
    peak = {}
    tot_peak = 0
    for rid in prof:
        tid = prof[rid]['tid']
        rprof = np.array(prof[rid]['rprofile'])
        if not valid_func(rprof): continue
        if len(rprof)-start-stop <= 0: continue
        peak_vec = identify_peaks(rprof, peak_func)
        # mask out regions outside of start and stop
        peak_vec[:start] = False
        if stop != 0:
            peak_vec[-stop:] = False
        peak[tid] = peak_vec
        tot_peak += sum(peak_vec)
    print "total valid transcrpts: {0} total peaks: {1}".format(len(peak), tot_peak)
    return peak

def get_peaks_from_histfile(fn, cds_range, start, stop):
    pb=parse_estimated_profile(fn)
    pc = codonp_from_basep(pb, cds_range, sum_frames)
    peak = get_peaks_from_profile(pc, start, stop, validate_profile, threshold)
    return peak

def unique_peaks(p1,p2):
    """ unique peaks in p1"""
    peak_unique = {}
    for tid, peak_vec1 in p1.iteritems():
        if tid not in p2:
            peak_unique[tid] = peak_vec1.copy()
        else: 
            peak_vec2 = p2[tid]
            peak_unique[tid] = peak_vec1 & (~peak_vec2)
    return peak_unique

def shared_peaks(p1,p2):
    """ shared peaks between p1 and p2"""
    return { tid: p1[tid] & p2[tid] for tid in set(p1.keys())&set(p2.keys()) }

def peak_sum(peak):
    return sum([sum(peak_vec) for _,peak_vec in peak.iteritems() ])

def single_peak_set_analysis(peak, start, stop, tseq, title, fn_prefix):
    codon_usage = get_codon_usage(peak.keys(), start, stop, tseq)
    codon_peak = get_codon_peak_cnt(peak, tseq)
    codon_freq = get_frequency(codon_peak, codon_usage)
    plot_codon_freq(codon_usage, title, fn_prefix+"_usage")
    plot_codon_freq(codon_peak, title, fn_prefix+"_peak")
    plot_codon_freq(codon_freq, title, fn_prefix+"_freq")
    aa_usage = codon_cnt_to_aa_cnt(codon_usage)
    aa_peak = codon_cnt_to_aa_cnt(codon_peak)
    aa_freq = get_frequency(aa_peak, aa_usage)
    plot_aa_freq(aa_usage, title, fn_prefix+"_usage")
    plot_aa_freq(aa_peak, title, fn_prefix+"_peak")
    plot_aa_freq(aa_freq, title, fn_prefix+"_freq")

def compare_two_peak_sets(peak1, peak2, start, stop, tseq, title, fn_prefix):
    single_peak_set_analysis(peak1, start, stop, tseq, title, fn_prefix+"1")
    single_peak_set_analysis(peak2, start, stop, tseq, title, fn_prefix+"2")
    p1_uniq = unique_peaks(peak1,peak2)
    p2_uniq = unique_peaks(peak2, peak1)
    p_share = shared_peaks(peak1, peak2)
    p1cnt = peak_sum(peak1)
    p2cnt = peak_sum(peak2)
    p1ucnt = peak_sum(p1_uniq)
    p2ucnt = peak_sum(p2_uniq)
    pscnt = peak_sum(p_share)
    single_peak_set_analysis(p1_uniq, start, stop, tseq, title, fn_prefix+"1uniq")    
    single_peak_set_analysis(p2_uniq, start, stop, tseq, title, fn_prefix+"2uniq")    
    single_peak_set_analysis(p_share, start, stop, tseq, title, fn_prefix+"_share")
    print "p1: unique:{0} ({1:.2%}) shared: {2} ({3:.2%}) total: {4}".format(p1ucnt, p1ucnt/float(p1cnt), pscnt, pscnt/float(p1cnt), p1cnt)
    print "p2: unique:{0} ({1:.2%}) shared: {2} ({3:.2%}) total: {4}".format(p2ucnt, p2ucnt/float(p2cnt), pscnt, pscnt/float(p2cnt), p2cnt)

def get_codon_usage(tid_list, start, stop, tseq):
    print "getting codon usage..."
    codon_cnt = {}
    for tid in tid_list:
        plen = len(tseq[tid])/3
        for i in xrange(start, plen-stop):
            codon = tseq[tid][i*3:(i+1)*3]
            codon_cnt.setdefault(codon,0)
            codon_cnt[codon]+=1
    return codon_cnt

def get_codon_peak_cnt(peak, tseq):
    print "getting codon peak count..."
    codon_cnt = {}
    for tid, p in peak.iteritems():
        for i in xrange(len(p)):
            if p[i] == True:
                codon = tseq[tid][i*3:(i+1)*3]
                codon_cnt.setdefault(codon,0)
                codon_cnt[codon] += 1
    return codon_cnt

def get_frequency(peak, usage):
    return { k: cnt/float(usage[k]) for k, cnt in peak.iteritems()}

def plot_codon_freq(codon_cnt, title, fn_prefix):
    codon_list = generate_codon_list()
    plot_cnt = len(codon_list)
    c = [ cmap(i) for i in np.linspace(0,1,plot_cnt) ]
    cnt_list = [ codon_cnt[codon] if codon in codon_cnt else 0 for codon in codon_list ]
    x = range(plot_cnt)
    plt.figure(figsize=(12,6))
    plt.bar(x, cnt_list, color=c, edgecolor='white',align='center')
    plt.xticks(x, codon_list, rotation='vertical', fontsize=12)
    plt.xlim((x[0]-1, x[-1]+1))
    plt.ylabel('peak count')
    plt.title(title)
    plt.savefig(fn_prefix+"_codon.pdf", bbox_inches='tight')
    plt.close()

def codon_cnt_to_aa_cnt(codon_cnt):
    aa_cnt = {}
    for codon, cnt in codon_cnt.iteritems():
        aa = codon2aa[codon]
        aa_cnt.setdefault(aa,0)
        aa_cnt[aa] += cnt
    return aa_cnt

def plot_aa_freq(aa_cnt, title, fn_prefix):
    aa_list = generate_aa_list()
    plot_cnt = len(aa_list)
    c = [ cmap(i) for i in np.linspace(0,1,plot_cnt) ]
    cnt_list = [ aa_cnt[aa] if aa in aa_cnt else 0 for aa in aa_list ]
    aa_name = [ aa2fullname[aa] for aa in aa_list]
    x = range(plot_cnt)
    plt.figure(figsize=(12,6))
    plt.bar(x, cnt_list, color=c, edgecolor='white',align='center')
    plt.xticks(x, aa_name, rotation='vertical', fontsize=12)
    plt.xlim((x[0]-1, x[-1]+1))
    plt.ylabel('peak count')
    plt.title(title)
    plt.savefig(fn_prefix+"_aa.pdf", bbox_inches='tight')
    plt.close()
    

if __name__ == "__main__": 
    cds_txt = "../ref/cds_range.txt"
    tfasta = "../ref/protein_coding_100_filtered.fasta"
    p_nc_nb_fn = "../ribomap/Lib-5-5-15_2_AGTTCC_R1_nodup.base"
    p_nc_wb_fn = "../ribomap/Lib-5-5-15_2_AGTTCC_R1_nonempty.base"
    p_wc_nb_fn = "../ribomap/Lib-5-5-15_6_TCCCGA_R1_nodup.base"
    p_wc_wb_fn = "../ribomap/Lib-5-5-15_6_TCCCGA_R1_nonempty.base"
    odir = "../figures/prelim/codon_usage/"
    start=20
    stop=20

    cds_range = get_cds_range(cds_txt)
    tseq = get_tseq(tfasta, cds_range)
    # no Chx collapse barcode
    p_ncnb = get_peaks_from_histfile(p_nc_nb_fn, cds_range, start, stop)
    # no Chx no collapse barcode
    p_ncwb = get_peaks_from_histfile(p_nc_wb_fn, cds_range, start, stop)
    # Chx collapse barcode
    p_wcnb = get_peaks_from_histfile(p_wc_nb_fn, cds_range, start, stop)
    # Chx no collapse barcode
    p_wcwb = get_peaks_from_histfile(p_wc_wb_fn, cds_range, start, stop)
    print "no chx barcode comparison 1: collapsed barcode 2:no collapse"
    compare_two_peak_sets(p_ncnb, p_ncwb, start, stop, tseq, 'no Chx', odir+'noChxYNbarcode')    
    print "10 chx barcode comarison 1: collapsed barcode 2: no collpase"
    compare_two_peak_sets(p_wcnb, p_wcwb, start, stop, tseq, '10x Chx', odir+'Chx10YNbarcode')        
    print "chx comarison (barcode collapsed) 1: no chx 2: 10 chx"
    compare_two_peak_sets(p_ncnb, p_wcnb, start, stop, tseq, 'barcode collapsed', odir+'BarcodeYNChx')        
    
    


