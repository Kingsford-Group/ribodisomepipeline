#!/usr/bin/env python
import numpy as np
import os
import sys
from footprint_hist_parser import *

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

def get_max_len(cds_range):
    max_len = 0
    for tid, (start, end) in cds_range.iteritems():
        if end > max_len:
            max_len = end
    return max_len

#=============================
# generate meta profile
#=============================
def offset_to_start(i, start):
    return i-start

def create_meta_profile(tlist, cds_range, rend, ibegin, iend):
    """ 
    create count accumulation for each base location with different read length
    end: whether 0-base is the START codon (5p) or STOP codon (3p)
    ibegin: index bound to include before START (negative)
    iend: bases to include after STOP (positive)
    """
    print "\ncreate meta hist"
    max_len = get_max_len(cds_range)
    meta_hist = {}
    if rend == "5p":
        iend = max_len + iend
    elif rend == "3p":
        ibegin = ibegin - max_len
    else:
        print "wrong rend parameter! can only take '5p' or '3p'"
        exit(1)
    for rid in tlist:
        tid = tlist[rid]['tid']
        start, end = cds_range[tid]
        if rend == "5p": anchor = start
        elif rend == "3p": anchor = end
        for rlen in tlist[rid]['prof']:
            for pos, cnt in tlist[rid]['prof'][rlen]:
                i = offset_to_start(pos, anchor)
                if i < ibegin or i > iend : continue
                meta_hist.setdefault(i,{})
                meta_hist[i].setdefault(rlen,0)
                meta_hist[i][rlen] += cnt
        sys.stdout.write("processed transcript {0}.\t\r".format(rid))
        sys.stdout.flush()
    sys.stdout.write("\n")
    return meta_hist

#===============================
# plot read length distribution
#===============================
def plot_read_len_hist(tlist, fn_prefix, rlen_min=0, rlen_max=1000):
    print "plotting read length distribution..."
    rlen2cnt = {}
    for rid in tlist:
        for rlen in tlist[rid]['prof']:
            tid = tlist[rid]['tid']
            for pos, cnt in tlist[rid]['prof'][rlen]:
                rlen2cnt.setdefault(rlen, 0)
                rlen2cnt[rlen] += cnt
    rlen_list = np.array(rlen2cnt.keys())
    select = (rlen_list >= rlen_min) & (rlen_list <= rlen_max)
    rlen_select = rlen_list[select]
    cnt_list = np.array(rlen2cnt.values())
    cnt_select = cnt_list[select]
    fig_width = len(rlen_list)/10.0*1.5
    plt.figure(figsize=(fig_width, 5))
    plt.bar(rlen_select-0.4, cnt_select, width=0.8, color='b', edgecolor='white', alpha=0.5)
    plt.xlabel('read length')
    plt.ylabel('read count')
    plt.xlim((min(rlen_select)-1, max(rlen_select)+1))
    plt.savefig(fn_prefix+"_rlen_hist.pdf", bbox_inches='tight')
    plt.close()

def plot_read_len_frame_hist(meta_hist, rlrange, rlen_shifts, ibegin, iend, norm, color, fn_prefix):
    width=int(50/10.0*1.5)
    plt.figure(figsize=(width, 6))
    print "plotting read length frame distribution..."
    rlmin, rlmax = rlrange
    rlmin -= 3
    rlmax += 3
    ymax = 0
    for i in meta_hist:
        if i<ibegin or i>iend : continue
        x_rlen = np.arange(rlmin, rlmax+1)
        y_cnt = np.array([0.1]*len(x_rlen))
        for rlen, cnt in meta_hist[i].iteritems():
            rlen += rlen_shifts[i%3]
            idx = rlen-rlmin
            y_cnt[idx] = cnt
        if norm == True:
            y_cnt /= np.mean(y_cnt[y_cnt!=0])
        y_cnt = np.log10(y_cnt)
        # plt.semilogy(rlen, cnt, color=color[i%3], lw=2, alpha=0.1)
        plt.plot(x_rlen, y_cnt, color=color[i%3], lw=2, alpha=0.1)
        ymax_loc = max(y_cnt)
        if ymax_loc > ymax:
            ymax = ymax_loc
    plt.vlines(range(1,rlmax,3), -1, ymax, colors='r', lw=2, alpha=0.3)
    plt.xlabel("read length")
    plt.ylabel("log footprint count")
    plt.xlim((rlmin,rlmax))
    plt.ylim((-1,ymax))
    xmin = rlmin - rlmin%5
    xmax = rlmax - rlmax%5
    plt.xticks(range(xmin, xmax+5, 5))
    # plt.show()
    plt.savefig("{0}_rlen_frame_hist.pdf".format(fn_prefix), bbox_inches='tight')
    plt.close()

#=============================
# main
#=============================
def plot_rlen_hist_pipe():
    if len(sys.argv) != 4:
        print "Usage: python read_len_hist.py cds_range_fname input_rlen.hist output_dir"
        exit(1)
    cds_txt = sys.argv[1]
    hist_fn = sys.argv[2]
    odir = sys.argv[3]
    utr5_offset = -24
    utr3_offset = 50
    imax = 300  # right most position form the 5' end to sample for historgram
    ensure_dir(odir)
    cds_range = get_cds_range(cds_txt)
    tlist = parse_rlen_hist(hist_fn)
    rlmin, rlmax = print_stats(tlist, cds_range)
    c = ['purple', 'blue', 'cyan' ]
    meta_hist = create_meta_profile(tlist, cds_range, "5p", utr5_offset, utr3_offset)
    fn_prefix = odir+"/"+get_file_core(hist_fn)
    plot_read_len_hist(tlist, fn_prefix, 33, 90)
    rlen_shifts = [ 0, 0, 0 ]
    # plot_rlen_hist(meta_hist, [rlmin, rlmax], rlen_shifts, utr5_offset, 0, False, c, fn_prefix+"_n{0}_ns".format(-utr5_offset))
    # plot_rlen_hist(meta_hist, [rlmin, rlmax], rlen_shifts, utr5_offset, imax, False, c, fn_prefix+"_p{0}_ns".format(imax))
    # rlen_shifts = [ 0, -2, -1]
    # plot_read_len_frame_hist(meta_hist, [rlmin, rlmax], rlen_shifts, 1, imax, False, c, fn_prefix)

if __name__ == "__main__": plot_rlen_hist_pipe()
