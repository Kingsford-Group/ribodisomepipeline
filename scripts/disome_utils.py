#!/usr/bin/env python
from footprint_hist_parser import parse_rlen_hist, get_cds_range

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

if __name__ == "__main__":
    cds_txt = "../ref/cds_range.txt"
    hist_fn = "../rlen_hist/prelim/Lib-5-5-15_3_CACGAT_R1_nodup.hist"
    cds_range = get_cds_range(cds_txt)
    tlist = parse_rlen_hist(hist_fn)
    
