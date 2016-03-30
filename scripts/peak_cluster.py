#!/usr/bin/env python
import numpy as np
from ribomap_result_parser import *

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['font.size'] = 20
rcParams['xtick.major.size'] = 5
rcParams['ytick.major.size'] = 5
cmap = matplotlib.cm.gist_rainbow

def validate_profile(vec, cnt_cutoff=1):
    # > 50% cnts > 1
    return np.mean(vec>cnt_cutoff)>0.5

def threshold(vec):
    return np.median(vec)+2*np.std(vec)

def threshold_with_min(vec, min_cnt=3):
    return max(min_cnt, threshold(vec))

def identify_peaks(vec, peak_func):
    return vec > peak_func(vec)

def cluster_peaks(peak, peak_width, seg_len):
    seg_list = []
    # round 1: find the first peak in the vector
    for i in xrange(len(peak)):
        if peak[i] == True:
            start = i
            end = i
            break
    # round 2: cluster peaks
    # cluster peaks within ribo_len
    i += 1
    while i < len(peak):
        if peak[i] == True:
            if i- end + 1 <= peak_width:
                end = i
            else:
                seg_list.append((start,end+1))
                start = i
                end = i
        i += 1
    # add last segment to list
    if len(seg_list)!=0 and end != seg_list[-1][1]:
        seg_list.append((start,end+1))
    return [ s for s in seg_list if (s[1]-s[0]+1)>seg_len ]

def peak_dist(p,d):
    """
    evidence that peaks tend to cluster together
    """
    pplen_list = []
    j = 0
    for tid in p:
        prof = p[tid]['rprofile']
        if not validate_profile(prof): continue
        peak = identify_peaks(prof, threshold)
        ploc = np.where(peak==True)[0]
        if len(ploc) < 2: continue
        pplen = [ ploc[i]-ploc[i-1] for i in xrange(1,len(ploc)) ]
        pplen_list.extend(pplen)
        j += 1
    print "number of transcripts included: {0}".format(j)
    print "average distance between peaks: {0}".format(np.mean(pplen_list))
    print "median distance between peaks: {0}".format(np.median(pplen_list))
    print "percentage within distance of {0}: {1:.2%}".format(d, np.mean(np.array(pplen_list)<d))
    plt.figure()
    ns,bins,patches = plt.hist(pplen_list,range(1,202,2))
    plt.xlabel("distance bwteen peaks")
    plt.ylabel("count")
    plt.savefig("peak_dist_hist.pdf", bbox_inches='tight')

def seg_len(seg_list):
    return [ s[1] - s[0] + 1 for s in seg_list ]

def jam_len(seg_list):
    return sum(seg_len(seg_list))

def jam_peak_cnt(peak, seg_list):
    return sum([np.sum(peak[s[0]:s[1]+1]) for s in seg_list])

def peak_segs(p,peak_width, seg_len):
    """
    stats on peak-rich regions
    """
    pcnt = 0
    jpcnt = 0
    clen = 0
    jclen = 0
    j = 0
    jratio_list = []
    for rid in p:
        prof = np.array(p[rid]['rprofile'])
        if not validate_profile(prof): continue
        peak = identify_peaks(prof, threshold)
        seg_list = cluster_peaks(peak,peak_width,seg_len)
        pcnt += np.sum(peak)
        clen += len(peak)
        this_jclen = jam_len(seg_list)
        jratio_list.append(float(this_jclen)/len(peak))
        jclen += this_jclen
        jpcnt += jam_peak_cnt(peak, seg_list)
        if j<10:
            # plot to eyeball bugs
            width=int(len(prof)/10.0*1.5)
            fig = plt.figure(figsize=(width,6))
            ax = fig.add_subplot(111)
            plt.plot(prof, 'b-', alpha=0.5)
            plt.scatter(np.where(peak==True)[0], prof[peak==True], color='r', marker='o')
            h = max(prof)*1.5
            for s in seg_list:
                plt.axvspan(s[0], s[1], color='r', alpha=0.2)
            plt.xlim((-1,len(prof)+2))
            plt.ylim((0,h))
            plt.title("{0}".format(p[rid]['tid']))
            #plt.show()
            plt.savefig("sample_profiles/{0}_peaks.pdf".format(p[rid]['tid']), bbox_inches='tight')
        # plt.close()
        j += 1
    print "number of transcripts included: {0}".format(j)
    print "total peaks: {0}".format(pcnt)
    print "clustered peaks: {0} ({1:.2%})".format(jpcnt, jpcnt/float(pcnt))
    print "total codon count: {0}".format(clen)
    print "peak-rich regions: {0} ({1:.2%})".format(jclen, jclen/float(clen))
    print "percentage of peaks in peak-rich regions: {0:.2%}".format(jpcnt/float(jclen))
    print "percentage of peaks {0:.2%}".format(pcnt/float(clen))
    plt.figure()
    ns,bins,patches = plt.hist(jratio_list,50)
    plt.xlabel("percentage of peak-rich regions")
    plt.ylabel("number of transcripts")
    plt.savefig("jam_portion.pdf", bbox_inches='tight')

def first_segs(p,peak_width, seg_len):
    """
    stats of  waterfall segs
    """
    j = 0
    slen = []
    epoint = []
    for tid in p:
        prof = p[tid]['rprofile']
        if not validate_profile(prof): continue
        peak = identify_peaks(prof, threshold)
        seg_list = cluster_peaks(peak,peak_width,seg_len)
        if len(seg_list)==0: continue
        if seg_list[0][0] > 50: continue
        slen.append(seg_list[0][1]-seg_list[0][0]+1)
        epoint.append(seg_list[0][1])
        j += 1
    print "number of transcripts included: {0}".format(j)    
    print "average waterfall segment length: {0}".format(np.mean(slen))
    plt.figure()
    ns,bins,patches = plt.hist(slen,20)
    plt.xlabel("waterfall segment length")
    plt.ylabel("number of transcripts")
    plt.savefig("first_seg_len_hist.pdf", bbox_inches='tight')
    print "average waterfall location: {0}".format(np.mean(epoint))
    plt.figure()
    ns,bins,patches = plt.hist(epoint, 20)
    plt.xlabel("waterfall location")
    plt.ylabel("number of transcripts")
    plt.savefig("first_seg_end_hist.pdf", bbox_inches='tight')

def main():
    p=parse_estimated_profile("BY_FP.codon")
    ribo_len = 10
    peak_segs(p, ribo_len*2, ribo_len)
    peak_dist(p, ribo_len*2)
    first_segs(p,ribo_len*2, ribo_len)

    j = 0
    for rid in p:
        tid = p[rid]['tid']
        rprof = np.array(p[rid]['rprofile'])
        if not validate_profile(rprof): continue
        peak = identify_peaks(rprof, threshold)
        seg_list = cluster_peaks(peak,ribo_len*2,ribo_len)
        j += 1
    
if __name__ == "__main__": main()
