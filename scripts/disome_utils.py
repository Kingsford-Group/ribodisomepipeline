#!/usr/bin/env python
from footprint_hist_parser import parse_rlen_hist, get_cds_range
from ribomap_result_parser import parse_estimated_profile

def generate_cds_coverage_profile(tlist, rlen_min, rlen_max, cds_range):
    p = {}
    for rid in tlist:
        tid = tlist[rid]['tid']
        start, stop = cds_range[tid]
        profile = np.zeros(stop-start)
        for rlen, pos_cnt_list in tlist[rid]['prof'].iteritems():
            if rlen < rlen_min or rlen > rlen_max: 
                continue
            else:
                for pos, cnt in pos_cnt_list:
                    for i in xrange(pos-start, pos-start+rlen):
                        if i<0 or i>= stop-start: continue
                        profile[i] += cnt
        if np.any(profile!=0):
            p[tid] = profile
    return p

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

def get_total_monosome_cnt(p, cds_range):
    tcnt = {}
    for rid in p:
        tid = p[rid]['tid']
        start, stop = cds_range[tid]
        tcnt[tid] = sum(p[rid]['rprofile'][start:stop])
    return tcnt

if __name__ == "__main__":
    cds_txt = "../ref/cds_range.txt"
    hist_fn = "../rlen_hist/prelim/Lib-5-5-15_3_CACGAT_R1_nodup.hist"
    sfn = "../ribomap/Lib-5-5-15_2_AGTTCC_R1_nodup.base"
    cds_range = get_cds_range(cds_txt)
    doublets = parse_rlen_hist(hist_fn)
    singlets = parse_estimated_profile(sfn)
    dcnt = get_total_disome_cnt(doublets, 40, 80, cds_range)
    scnt = get_total_monosome_cnt(singlets, cds_range)
    collision_ratio = []
    for tid in scnt:
        if tid in dcnt:
            ratio = float(dcnt[tid])/(dcnt[tid] + scnt[tid])
        else:
            ratio = 0
        collision_ratio.append(ratio)
    
    import matplotlib.pyplot as plt
    import numpy as np
    ns,bins,patches = plt.hist(collision_ratio, np.linspace(0,1,50), log=True)
    plt.savefig('test.pdf', bbox_inches='tight')
    plt.close()
