#!/usr/bin/env python
import numpy as np
from doublet_profile import generate_codon_profile_from_rlen_hist, get_tid2codonp_from_ribomap_base
from significant_doublet_count import get_window_cnt, read_pvals_from_file
from footprint_hist_parser import get_cds_range, get_tseq
from file_names import *

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['font.size'] = 12
rcParams['xtick.major.size'] = 5
rcParams['ytick.major.size'] = 5

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

def get_values(tid2pos, profile_list, func=None):
    vec = []
    for tid in tid2pos:
        for pos in tid2pos[tid]:
            if func == None:
                vec.append(profile_list[tid][pos])
            else:
                vec.append(func(pos, profile_list[tid]))
    return vec

def get_collision_rates(tid2pos, sprof, dprof):
    vec = []
    for tid in tid2pos:
        for pos in tid2pos[tid]:
            dcnt = get_window_cnt(pos, dprof[tid])
            scnt = sprof[tid][pos]
            cr = float(dcnt)/scnt
            vec.append(cr)
    return vec

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
    tid2pvals = read_pvals_from_file(pfname)
    tid2sig, tid2rest = get_significant_peaks(tid2pvals)
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

if __name__ == "__main__":
    cds_range = get_cds_range(cds_txt) 
    significance_analysis_pipeline(nchx_pfn, nchx_dfn, nchx_sfn, "../ds_cmp/nchx_pval")
    significance_analysis_pipeline(chx_pfn, chx_dfn, chx_sfn, "../ds_cmp/chx_pval")
    significance_analysis_pipeline(wt_pfn, wt_dfn, wt_sfn, "../ds_cmp/wt_pval")
    significance_analysis_pipeline(dom34_pfn, dom34_dfn, dom34_sfn, "../ds_cmp/dom34_pval")

