#!/usr/bin/env python
import numpy as np
import os
import sys
from io_utils import *

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['font.size'] = 12
rcParams['xtick.major.size'] = 5
rcParams['ytick.major.size'] = 5

#=============================
# utils
#=============================
def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)

def get_file_core(fname):
    istart = fname.rfind("/")
    iend = fname.rfind(".hist")
    return fname[istart:iend]

def unzip_dic_tuples(dic):
    keys, vals = zip(*dic.items())
    order = np.argsort(keys)
    return np.array(keys)[order], np.array(vals)[order]

def get_max_len(cds_range):
    max_len = 0
    for tid, (start, end) in cds_range.iteritems():
        if end > max_len:
            max_len = end
    return max_len
    
def is_sparse(vec, threshold):
    order = np.argsort(vec)
    # largest peak not included when determining sparsity
    if np.sum(vec[order][:-1])>= threshold:
        return False
    else:
        return True

def get_ylim(vec, pcnt=3):
    """ 
    pcnt: which peak to use as the ylim
    """
    order = np.argsort(vec)
    idx = order[-pcnt]
    return vec[idx]

def get_pos_hist(pos_hist, start, end):
    y_cnt = np.array([0]*(end-start+1))
    for pos, cnt in pos_hist.iteritems():
        if pos < start or pos > end: continue
        idx = pos - start
        y_cnt[idx] = cnt
    return y_cnt

def get_pobs(pos_hist, rlstart, rlend, istart, iend):
    pobs = {}
    for rlen in xrange(rlstart, rlend+1):
        pobs[rlen] = get_pos_hist(pos_hist[rlen], istart, iend)
    return pobs

def get_frames(vec):
    fsum = np.array([ np.sum(vec[i::3]) for i in xrange(3) ])
    tsum = float(sum(fsum))
    f0 = np.argmax(fsum)
    fdic = { f0:0, (f0+1)%3:1, (f0+2)%3:2 }
    print "frames: {0:.2%} {1:.2%} {2:.2%}".format(fsum[0]/tsum, fsum[1]/tsum, fsum[2]/tsum)
    return fdic

def get_frame_str_vec(vec):
    fsum = np.array([ np.sum(vec[i::3]) for i in xrange(3) ])
    tsum = float(sum(fsum))
    f0 = np.argmax(fsum)
    fdic = { f0:0, (f0+1)%3:1, (f0+2)%3:2 }
    sframe =  "{0:.2%} {1:.2%} {2:.2%}".format(fsum[0]/tsum, fsum[1]/tsum, fsum[2]/tsum)
    return sframe  

def get_frame_str(pobs):
    return { rlen: get_frame_str_vec(vec) for rlen, vec in pobs.iteritems() }

#=============================
# transcript pre-selection
#=============================
def filter_transcript_by_length(cds_range, length_cutoff):
    tid_list = np.array(cds_range.keys())
    tlen = np.array([ v[1]-v[0] for v in cds_range.values()])
    return set(tid_list[tlen > length_cutoff])

#=============================
# generate meta profile
#=============================
def create_rlen_meta_profile(tlist, cds_range, tid_included, ibegin, iend):
    """ 
    create count accumulation for each base location with different read length
    exclude transcript if id not in tid_included (a set)
    ibegin: index bound to include before START (negative)
    iend: number of bases to include after START (positive)
    """
    print "\ncreate meta hist"
    meta_hist = {}
    for rid in tlist:
        tid = tlist[rid]['tid']
        if tid not in tid_included: continue
        start, end = cds_range[tid]
        for rlen in tlist[rid]['prof']:
            meta_hist.setdefault(rlen, {})
            for pos, cnt in tlist[rid]['prof'][rlen]:
                i = pos - start
                if i < ibegin or i > iend : continue
                meta_hist[rlen].setdefault(i,0)
                meta_hist[rlen][i] += cnt
        sys.stdout.write("processed transcript {0}.\t\r".format(rid))
        sys.stdout.flush()
    sys.stdout.write("\n")
    return meta_hist

#=============================
# plot meta profile
#=============================
def plot_pos_hist(meta_hist, ibegin, iend, min_cnt, fn_prefix):
    print "plotting pos hist..."
    rlen_list = sorted(meta_hist.keys())
    x_pos = np.arange(ibegin, iend+1)
    figwidth = len(x_pos)/10.0*1.5
    for rlen in rlen_list:
        y_cnt = get_pos_hist(meta_hist[rlen], ibegin, iend)
        if is_sparse(y_cnt, min_cnt): continue
        fig = plt.figure(figsize=(figwidth,6))
        ax = fig.add_subplot(1,1,1)
        max_cnt = get_ylim(y_cnt, pcnt=1)
        ymax = max_cnt + 1
        plt.bar(x_pos-0.4, y_cnt, width=0.8, color='b', edgecolor='white', alpha=0.5)
        # plt.bar(-rlen+16-0.4, ymax, width=0.8, color='r', edgecolor='white', alpha=0.3)
        plt.xlabel("offset from start")
        plt.ylabel("footprint count")
        plt.xticks(range(-100, iend+1, 5))
        plt.xlim((ibegin-1, iend+1))
        plt.ylim((0.1, ymax))
        frame_str = get_frame_str_vec(y_cnt)
        plt.title("read length = {0}\nframe distribution: {1}".format(rlen, frame_str), size=20)
        plt.savefig("{0}_{1}_pos_hist.pdf".format(fn_prefix, rlen), bbox_inches="tight")
        plt.close()

def plot_meta_pos_hist(meta_hist, ibegin, iend, fn_prefix):
    print "plotting meta pos hist..."
    x_pos = np.arange(ibegin, iend+1)
    y_cnt = np.array([0]*len(x_pos))
    figwidth = len(x_pos)/10.0*1.5
    for rlen in meta_hist:
        y_cnt_rlen = get_pos_hist(meta_hist[rlen], ibegin, iend)
        y_cnt += y_cnt_rlen
    fig = plt.figure(figsize=(figwidth,6))        
    ax = fig.add_subplot(1,1,1)
    max_cnt = get_ylim(y_cnt)
    ymax = max_cnt + 1
    plt.bar(x_pos-0.4, y_cnt, width=0.8, color='b', edgecolor='white', alpha=0.5)
    imax = np.argmax(y_cnt)
    plt.xlabel("offset from start")
    plt.ylabel("footprint count")
    plt.xticks(range(-100, iend+1, 5))
    plt.xlim((ibegin-1, iend+1))
    plt.ylim((0.1, ymax))
    plt.savefig("{0}_{1}_pos_hist.pdf".format(fn_prefix, "meta"), bbox_inches="tight")
    plt.close()

#=============================
# main
#=============================
def plot_rlen_hist_pipe():
    if len(sys.argv) != 4:
        print "Usage: python meta_profile.py cds_range_fname input_rlen.hist output_dir"
        exit(1)
    cds_txt = sys.argv[1]
    hist_fn = sys.argv[2]
    odir = sys.argv[3]
    utr5_offset = -24
    imax = 300  # right most position form the 5' end to sample for historgram
    min_sample_cnt = 100
    ensure_dir(odir)
    cds_range = get_cds_range(cds_txt)
    tid_select = filter_transcript_by_length(cds_range, imax)
    tlist = parse_rlen_hist(hist_fn)
    rlmin, rlmax = print_stats(tlist, cds_range)
    meta_hist = create_rlen_meta_profile(tlist, cds_range, tid_select, utr5_offset, imax)
    fn_prefix = odir+"/"+get_file_core(hist_fn)
    plot_pos_hist(meta_hist, utr5_offset, imax, min_sample_cnt, fn_prefix)    
    plot_meta_pos_hist(meta_hist, utr5_offset, imax, fn_prefix)

if __name__ == "__main__": plot_rlen_hist_pipe()
