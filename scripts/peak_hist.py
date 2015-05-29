#!/usr/bin/env python
import sys
import numpy as np
import scipy.stats
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

#def main():
def get_codon_usage(tid_list, start, stop, tseq):
    print "getting codon usage..."
    codon_cnt = {}
    for tid in tid_list:
        plen = len(tseq[tid])/3
        for i in xrange(start, plen-stop):
            codon = tseq[tid][i:i+3]
            codon_cnt.setdefault(codon,0)
            codon_cnt[codon]+=1
    return codon_cnt

def get_codon_peak_cnt(peak, tseq):
    print "getting codon peak count..."
    codon_cnt = {}
    for tid, p in peak.iteritems():
        for i in xrange(len(p)):
            if p[i] == True:
                codon = tseq[tid][i:i+3]
                codon_cnt.setdefault(codon,0)
                codon_cnt[codon] += 1
    return codon_cnt

def get_codon_frequency(codon_peak, codon_usage):
    return { codon: cnt/float(codon_usage[codon]) for codon, cnt in codon_peak.iteritems()}

def generate_codon_list():
    codon_list = []
    base_list = ['G', 'A', 'C', 'T']
    stop_list = ['TAA', 'TAG', 'TGA']
    for b0 in base_list:
        for b1 in base_list:
            for b2 in base_list:
                codon = b0+b1+b2
                if codon not in stop_list:
                    codon_list.append(codon)
    return codon_list

def plot_codon_freq(codon_cnt, title, fn):
    codon_list = generate_codon_list()
    plot_cnt = len(codon_list)
    c = [ cmap(i) for i in np.linspace(0,1,plot_cnt) ]
    cnt_list = [ codon_cnt[codon] for codon in codon_list]
    x = range(plot_cnt)
    plt.figure(figsize=(12,6))
    plt.bar(x, cnt_list, color=c, edgecolor='white',align='center')
    plt.xticks(x, codon_list, rotation='vertical', fontsize=12)
    plt.xlim((x[0]-1, x[-1]+1))
    plt.savefig(fn, bbox_inches='tight')
    plt.close()


if __name__ == "__main__": 
    if len(sys.argv)!=4:
        print "Usage: python peak_hist.py cds_range.txt ref.fasta ribomap.base"
        exit(1)

    cds_txt = sys.argv[1]
    tfasta = sys.argv[2]
    pfn = sys.argv[3]
    cds_range = get_cds_range(cds_txt)
    tseq = get_tseq(tfasta, cds_range)
    pb=parse_estimated_profile(pfn)
    pc = codonp_from_basep(pb, cds_range, sum_frames)
    rc_threshold = 1
    start=20
    stop=20
    peak = get_peaks_from_profile(pc, start, stop, validate_profile, threshold)
    codon_usage = get_codon_usage(peak.keys(), start, stop, tseq)
    codon_peak = get_codon_peak_cnt(peak, tseq)
    codon_freq = get_codon_frequency(codon_peak, codon_usage)
    plot_codon_freq(codon_freq, 'test', 'test.pdf')
    


