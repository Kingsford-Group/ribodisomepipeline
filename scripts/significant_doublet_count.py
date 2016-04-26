#!/usr/bin/env python
import os
import sys
import numpy as np
from doublet_profile import *

def single_peak_pval(peak_pos, doublet_profile, sd_distance=14, window_size=10, sample_size=10000):
    """ 
    empirical p-value calculation: probability of observing doublets more than expected by chance
    sd_distance: upstream distance between doublet and singlet (unit: codon)
    window_size: total number of positions for collecting doublet count (unit: codon)
    """
    iend = peak_pos - sd_distance
    istart = iend - window_size
    if istart < 0: istart = 0
    npos = iend - istart 
    dcnt_fg = np.sum(doublet_profile[istart:iend])
    dcnt_bg = np.array([ np.sum(np.random.choice(doublet_profile, npos, replace=False)) for i in xrange(sample_size) ])
    return np.mean(dcnt_fg <= dcnt_bg)

def get_window_cnt(peak_pos, profile, sd_distance=sd_distance, window_size=window_size):
    iend = peak_pos - sd_distance
    if iend < 0 : return None
    istart = iend - window_size
    if istart >= len(profile): return None
    if istart < 0: istart = 0
    return np.sum(profile[istart:iend])

def batch_pval_per_transcript(peak_vec, profile, sd_distance=sd_distance, window_size=window_size, sample_size=sample_size):
    pval_list = []
    dcnt_bg = np.array([ np.sum(np.random.choice(profile, window_size, replace=False)) for i in xrange(sample_size) ])
    for peak_pos in np.where(peak_vec)[0]:
        dcnt_fg = get_window_cnt(peak_pos, profile, sd_distance, window_size)
        if dcnt_fg == None: continue
        pval = np.mean(dcnt_fg <= dcnt_bg)
        pval_list.append((peak_pos, pval))
    return pval_list

def batch_peak_pval(tid2peaks, profiles, sd_distance, window_size):
    tid2pvals = {}
    for tid in tid2peaks:
        if tid not in profiles: continue
        pval_list = batch_pval_per_transcript(tid2peaks[tid], profiles[tid], sd_distance, window_size)
        if len(pval_list) != 0: tid2pvals[tid] = pval_list
    return tid2pvals

def write_pvals_to_file(tid2pvals, fname):
    tf = open(fname, 'w')
    for tid in tid2pvals:
        words = [ "{0}:".format(tid) ]
        for pos, pval in tid2pvals[tid]:
            words.append( "{0}-{1}".format(pos,pval))
        tf.write("\t".join(words)+"\n")
    tf.close()

def read_pvals_from_file(fname):
    tid2pvals = {}
    tf = open(fname)
    for line in tf:
        words = line.rstrip('\n').split('\t')
        tid = words[0].rstrip(':')
        pval_list = []
        for spair in words[1:]:
            sep = spair.find('-')
            pval_list.append((int(spair[:sep]), float(spair[sep+1:])))
        tid2pvals[tid] = pval_list
    tf.close()
    return tid2pvals

def collision_significance_pipe(cds_range, sfname, dfname, ofname, sd_distance, window_size, peak_type='singlet'):
    print "doublet", os.path.basename(dfname)
    dcp = generate_codon_profile_from_rlen_hist(dfname, cds_range)
    print "singlet", os.path.basename(sfname)
    scp = get_tid2codonp_from_ribomap_base(sfname, cds_range)
    tid_list = get_hc_list_from_prof(scp)
    if peak_type == 'singlet':
        profiles = dcp
    elif peak_type == 'doublet':
        profiles = scp
    else:
        print "peak type {0} not supported!".format(peak_type)
        exit(1)
    peaks = batch_peak_call(scp)
    peak_hc = { tid:peaks[tid] for tid in tid_list if tid in peaks }
    tid2pvals = batch_peak_pval(peak_hc, profiles, sd_distance, window_size)
    write_pvals_to_file(tid2pvals, ofname)

if __name__ == "__main__":
    cds_range = get_cds_range(cds_txt)
    if len(sys.argv) != 7:
        print "Usage: python significant_doublet_count.py doublet.hist singlet.ribomap output.txt"
        exit(1)
    dfname = sys.argv[1]
    sfname = sys.argv[2]
    ofname = sys.argv[3]
    sd_distance = int(sys.argv[4])
    window_size = int(sys.argv[5])
    peak_type = sys.argv[6]

    collision_significance_pipe(cds_range, sfname, dfname, ofname, sd_distance, window_size, peak_type)    
