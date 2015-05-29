#!/usr/bin/env python
from ribomap_result_parser import *
import sys

def frame_stats(p,cds_range, afterstart, beforestop):
    frame = [0]*3
    for rid in p:
        tid = p[rid]['tid']
        start,stop = cds_range[tid]
        start += afterstart
        if start > stop: continue
        stop -= beforestop
        if stop < start: continue
        rprof = p[rid]['rprofile'][start:stop]
        for i in xrange(3):
            frame[i] += sum(rprof[i::3])
    fsum = float(sum(frame))
    print "frame cnts: {0:.2f} {1:.2f} {2:.2f} total: {3:.2f}".format(frame[0], frame[1], frame[2], fsum)
    print "frames portion: {0:.2%} {1:.2%} {2:.2%}".format(frame[0]/fsum, frame[1]/fsum, frame[2]/fsum)
    

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print "Usage: python mapping_stats.py cds_range.txt ribomap.base"
        exit(1)
    
    afterstart=0 #20*6
    beforestop=0 #20*6
    cds_txt = sys.argv[1]
    pfn = sys.argv[2]
    cds_range = get_cds_range(cds_txt)
    p=parse_estimated_profile(pfn)
    frame_stats(p,cds_range, afterstart, beforestop)
    
