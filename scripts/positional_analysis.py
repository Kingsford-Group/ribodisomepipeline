#!/usr/bin/env python
import numpy as np
from doublet_profile import generate_codon_profile_from_rlen_hist, get_tid2codonp_from_ribomap_base, get_codon_profile_from_deblur_profile
from io_utils import get_cds_range
from file_names import *

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['font.size'] = 12
rcParams['xtick.major.size'] = 5
rcParams['ytick.major.size'] = 5

def get_accumulative_count_relative_pos(tid2prof, num_points, len_lim=200, normed=True):
    bin_cnt = np.zeros(num_points)
    tcnt = 0
    for tid, prof in tid2prof.iteritems():
        # if len(prof)< len_lim: continue
        tcnt += 1
        interval = len(prof)/float(num_points)
        for i in xrange(len(prof)):
            bi = int(np.floor(i/interval))
            bin_cnt[bi] += prof[i]
    print "total read count:{0:.0f}".format(sum(bin_cnt))
    if normed:
        bin_cnt /= sum(bin_cnt)
    acc_cnt = np.zeros(num_points+1)
    for i in xrange(num_points):
        acc_cnt[i+1] = acc_cnt[i] + bin_cnt[i]
    # print "included transcripts {0} ({1:.2%})".format(tcnt, float(tcnt)/len(tid2prof))
    return acc_cnt

def get_longest_transcript(tid2prof):
    return max(map(len, tid2prof.values()))

def get_accumulative_count_absolute_pos(tid2prof, len_lim=200, normed=True):
    # initialize count vector
    if len_lim == None:
        num_points = get_longest_transcript(tid2prof)
    else:
        num_points = len_lim
    bin_cnt = np.zeros(num_points)
    # accumulate counts
    tcnt = 0
    for tid, prof in tid2prof.iteritems():
        # if len(prof) < num_points: continue
        tcnt += 1
        n = min(len(prof), num_points)
        bin_cnt[:n] += prof[:n]
    print "length included:{0} total read count:{1:.0f}".format(num_points, sum(bin_cnt))
    # normalize counts
    if normed:
        bin_cnt /= sum(bin_cnt)
    # construct cdf
    acc_cnt = np.zeros(num_points+1)
    for i in xrange(num_points):
        acc_cnt[i+1] = acc_cnt[i] + bin_cnt[i]
    # print "included transcripts {0} ({1:.2%})".format(tcnt, float(tcnt)/len(tid2prof))
    return acc_cnt

def positional_analysis_pipeline(sfname, dfname, cds_range, len_lim, oprfx, multimap=False):
    dcp = generate_codon_profile_from_rlen_hist(dfname, cds_range)
    if multimap == True:
        scp = get_tid2codonp_from_ribomap_base(sfname, cds_range)
    else:
        scp = get_codon_profile_from_deblur_profile(sfname)
    accd_rpos = get_accumulative_count_absolute_pos(dcp, len_lim)
    accs_rpos = get_accumulative_count_absolute_pos(scp, len_lim)

    plt.figure()
    plt.plot(accd_rpos, c='b', lw=1.5)
    plt.plot(accs_rpos, c='r', lw=1.5, ls='--')
    plt.legend(['doublet', 'singlet'], loc='upper left', frameon=False)
    plt.title('absolute positional cdf')
    plt.xlabel('position')
    plt.ylabel('frequency of included ribosome footprint')
    plt.ylim((0,1))
    plt.savefig('{0}_absolute_pos_cdf_{1}.pdf'.format(oprfx, len_lim))
    plt.close()

    accd_rpos = get_accumulative_count_relative_pos(dcp, 100, len_lim)
    accs_rpos = get_accumulative_count_relative_pos(scp, 100, len_lim)
    plt.figure()
    plt.plot(accd_rpos, c='b', lw=1.5)
    plt.plot(accs_rpos, c='r', lw=1.5, ls='--')
    plt.legend(['doublet', 'singlet'], loc='upper left', frameon=False)
    plt.title('relative positional cdf')
    plt.xlabel('relative position')
    plt.ylabel('frequency of included ribosome footprint')
    plt.ylim((0,1))
    plt.savefig('{0}_relative_pos_cdf_{1}.pdf'.format(oprfx, len_lim))
    plt.close()

if __name__ == "__main__":
    cds_range = get_cds_range(cds_txt)
    len_lim = 200
    positional_analysis_pipeline(nchx_sfn, nchx_dfn, cds_range, len_lim, "../uniquely_mapped/figures/nchx", multimap)
    positional_analysis_pipeline(chx_sfn, chx_dfn, cds_range, len_lim, "../uniquely_mapped/figures/chx", multimap)
    positional_analysis_pipeline(wt_sfn, wt_dfn, cds_range, len_lim, "../uniquely_mapped/figures/wt", multimap)
    positional_analysis_pipeline(dom34_sfn, dom34_dfn, cds_range, len_lim, "../uniquely_mapped/figures/dom34", multimap)
