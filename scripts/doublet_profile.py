#!/usr/bin/env python
import os
import sys
import numpy as np
import scipy.stats
from footprint_hist_parser import parse_rlen_hist, get_cds_range, get_tseq
from ribomap_result_parser import parse_estimated_profile
from ribofit_utils import validate_profile, codonp_from_basep, sum_frames
from peak_cluster import threshold, threshold_with_min, identify_peaks, cluster_peaks
from file_names import *

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['font.size'] = 12
rcParams['xtick.major.size'] = 5
rcParams['ytick.major.size'] = 5

def generate_cds_coverage_profile(tlist, rlen_min, rlen_max, cds_range):
    """ 
    coverage profile within cds for with length [rlen_min, rlen_max]
    """
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

def generate_cds_profile(tlist, rlen_min, rlen_max, cds_range, offset=0):
    """
    read start pileups within cds with length [rlen_min, rlen_max]
    """
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
                    i = pos + offset
                    if i<0 or i>= stop-start: continue
                    profile[i] += cnt
        if np.any(profile!=0):
            p[tid] = profile
    return p

def get_tid2basep_from_ribomap_base(ribomap_fname, cds_range):
    """ dict of tid : nucleotide profile (cds only) """
    sp_base = parse_estimated_profile(ribomap_fname)
    sp_cds = {}
    for rid in sp_base:
        tid = sp_base[rid]['tid']
        start, stop = cds_range[tid]
        pcds = sp_base[rid]['rprofile'][start:stop]
        sp_cds[tid] = np.array(pcds)
    return sp_cds

def merge_count_per_codon(tid2prof, merge_func):
    """
    input: dic tid : nucleotide profile
    output: dic tid : codon profile
    """
    return { tid : np.array(merge_func(p)) for tid, p in tid2prof.iteritems() }

def batch_peak_call(tid2prof, peak_min_cutoff=3):
    thresh = lambda vec : threshold_with_min(vec, peak_min_cutoff)
    tid2peak = { tid: identify_peaks(prof, thresh) for tid, prof in tid2prof.iteritems() }
    return { tid : peak for tid, peak in tid2peak.iteritems() if np.any(peak!=False) }
        
def get_tid2codonp_from_ribomap_base(ribomap_fname, cds_range, merge_func=sum_frames):
    """
    input: ribomap nucleotide profile
    output: dic -- transcript id : codon profile
    """
    # singlets (dic of ribomap outputs: rid: tid, rprofile, nprofile, mprofile)
    sp_base = parse_estimated_profile(ribomap_fname)
    sp_codon = codonp_from_basep(sp_base, cds_range, merge_func)
    return { p['tid'] : p['rprofile'] for rid, p in sp_codon.iteritems() }

def get_hc_list_from_prof(prof):
    """ get high coverage transcript list """
    return [ tid for tid in prof if validate_profile(prof[tid]) ]

def get_peak_segs(codonp, peak_width=3, peak_min_cutoff=3):
    """ output dic: transcript id : list of peak locations (start, end, count) """
    print "calling peaks..."
    thresh = lambda vec : threshold_with_min(vec, peak_min_cutoff)
    peak_segs = {}
    tcnt = 0
    pcnt = 0
    for tid, prof in codonp.iteritems():
        peak = identify_peaks(prof, thresh)
        if np.all(peak==False): continue
        seg_list = cluster_peaks(peak, peak_width, 0)
        if len(seg_list)==0: continue
        seg_cnt_list = [ (start,end,np.sum(prof[start:end])) for start,end in seg_list]
        peak_segs[tid] = seg_cnt_list
        tcnt += 1
        pcnt += len(seg_list)
        sys.stdout.write("processed transcript {0}.\t\r".format(tcnt))
        sys.stdout.flush()
    sys.stdout.write('\n')
    print "total peaks (merged) {0}".format(pcnt)
    return peak_segs

def generate_base_profile_from_rlen_hist(fname, cds_range):
    """
    default globals: dlen_min, dlen_max, peak_width, peak_min_cutoff
    """
    print os.path.basename(fname)
    reads = parse_rlen_hist(fname)
    basep = generate_cds_profile(reads, dlen_min, dlen_max, cds_range)
    return basep

def generate_codon_profile_from_rlen_hist(fname, cds_range):
    basep = generate_base_profile_from_rlen_hist(fname, cds_range)
    codonp = merge_count_per_codon(basep, sum_frames)
    return codonp

def peak_segs_from_rlen_hist_pipeline(fname, cds_range):
    """
    build a list of peak start,end,count per transcript from rlen hist
    """
    codonp = generate_codon_profile_from_rlen_hist(fname, cds_range)
    peak_segs = get_peak_segs(codonp, peak_width, peak_min_cutoff)
    return peak_segs

def get_top_tid_list(prof, top_func=np.mean, percentile=0.1):
    tid_list = np.array(prof.keys())
    dd = np.array([ top_func(prof[tid]) for tid in tid_list ])
    order = np.argsort(dd)[::-1]
    percent_cutoff = int(len(order)*percentile)
    return tid_list[order][:percent_cutoff].flatten()

def plot_dsprof(dp, sp, tid, figname, grid=False):
    cds_len = len(dp)
    assert len(dp) == len(sp)
    x_pos = np.arange(cds_len)
    xgrid = np.arange(0,cds_len, 3)
    width = min(300, int(cds_len/10.0*1.5))
    plt.figure(figsize=(width, 15))
    plt.subplot(2,1,1)
    if grid == True:
        plt.vlines(xgrid, [0]*len(xgrid), [max(dp)+1]*len(xgrid), 'r', linewidths=5, alpha=0.1)
    plt.vlines(x_pos, [0]*cds_len, dp, 'k', linewidths=5, alpha=0.5)
    plt.xticks(np.arange(0,cds_len, 10))
    plt.xlim((-1, cds_len))
    plt.ylim((0,max(dp)+1))
    plt.title("{0}\ndoublet profile".format(tid))
    plt.subplot(2,1,2)
    if grid == True:
        plt.vlines(xgrid, [0]*len(xgrid), [max(sp)+1]*len(xgrid), 'r', linewidths=5, alpha=0.1)
    plt.vlines(x_pos, [0]*cds_len, sp, 'k', linewidths=5, alpha=0.5)
    plt.xticks(np.arange(0,cds_len, 10))
    plt.xlim((-1, cds_len))
    plt.ylim((0, max(sp)+1))
    cr = float(sum(dp))/(sum(dp)+sum(sp))
    plt.title("collision rate {0:.2%}\nsinglet profile".format(cr))
    plt.xlabel("position")
    try:
        plt.savefig(figname, bbox_inches='tight')
    except ValueError:
        print "render error occurs! try without tight layout"
        try: 
            plt.savefig(figname, bbox_inches='tight')
        except ValueError:
            print "render error occurs without tight layout!"
    plt.close()
    print "{0} len: {1} average doublet: {2:.2f} collision rate: {3:.2%}".format(tid, len(sp), np.mean(dp), cr)

def plot_high_coverage_doublets(cds_range, fig_dir):
    dpb_nchx = generate_base_profile_from_rlen_hist(nchx_dfn, cds_range)
    dpb_chx = generate_base_profile_from_rlen_hist(chx_dfn, cds_range)
    dpb_wt = generate_base_profile_from_rlen_hist(wt_dfn, cds_range)
    dpb_dom34 = generate_base_profile_from_rlen_hist(dom34_dfn, cds_range)

    tid_dset = set(get_top_tid_list(dpb_nchx, np.mean, 0.01))
    tid_dset &= set(get_top_tid_list(dpb_chx, np.mean, 0.01))
    tid_dset &= set(get_top_tid_list(dpb_wt, np.mean, 0.01))
    tid_dset &= set(get_top_tid_list(dpb_dom34, np.mean, 0.01))
    print "# high coverage doublet profiles {0}".format(len(tid_dset))

    spb_nchx = get_tid2basep_from_ribomap_base(nchx_sfn, cds_range)
    spb_chx = get_tid2basep_from_ribomap_base(chx_sfn, cds_range)
    spb_wt = get_tid2basep_from_ribomap_base(wt_sfn, cds_range)
    spb_dom34 = get_tid2basep_from_ribomap_base(dom34_sfn, cds_range)

    tid_sset = set(spb_nchx.keys())
    tid_sset &= set(spb_chx.keys())
    tid_sset &= set(spb_wt.keys())
    tid_sset &= set(spb_dom34.keys())

    tid_share = list(tid_dset & tid_sset)
    print len(tid_share)
    for tid in tid_share:
        plot_dsprof(dpb_chx[tid], spb_chx[tid], "{0}_chx".format(tid), "{0}{1}_1chx.pdf".format(fig_dir,tid))
        plot_dsprof(dpb_nchx[tid], spb_nchx[tid], "{0}_nchx".format(tid), "{0}{1}_2nchx.pdf".format(fig_dir,tid))
        plot_dsprof(dpb_wt[tid], spb_wt[tid], "{0}_wt".format(tid), "{0}{1}_3wt.pdf".format(fig_dir,tid))
        plot_dsprof(dpb_dom34[tid], spb_dom34[tid], "{0}_dom34".format(tid), "{0}{1}_4dom34.pdf".format(fig_dir,tid))


def generate_dist_list(dpeaks, speaks, dstart, dend):
    dist_list = []
    dcnt = 0
    for tid in dpeaks:
        dcnt += np.sum(dpeaks[tid])
        if tid not in speaks: continue
        for dpos in np.where(dpeaks[tid])[0]:
            istart = dpos+dstart
            iend = dpos+dend
            if istart<0 or istart >= len(dpeaks[tid]): continue
            if dstart > 0 and dend > dstart:
                spos_list = np.where(speaks[tid][istart:iend])[0]
                if len(spos_list)== 0: continue
                spos = spos_list[0] + istart 
            if dstart <= 0 and dend < dstart:
                spos_list = np.where(speaks[tid][iend:istart])[0]
                if len(spos_list) == 0: continue
                spos = spos_list[-1]
            if (dstart>0 and dend<dstart) or (dstart<=0 and dend>dstart):
                print "wrong definition of dstart {0} and dend {1}!".format(dstart, dend)
                exit(1)
            dist_list.append([tid, dpos, spos])
    print "coupled peaks: {0} out of {1} ({2:.2%})".format(len(dist_list), dcnt, len(dist_list)/float(dcnt))
    return dist_list

def get_doublet_singlet_peak_base_dist(dfname, sfname, cds_range, dstart=40, dend=100):
    """
    return list: [transcript_id, doublet_loc, doublet_cnt, singlet_loc, singlet_cnt]
    default globals: dlen_min, dlen_max, peak_width, peak_min_cutoff
    """
    dbp = generate_base_profile_from_rlen_hist(dfname, cds_range)
    dpeaks = batch_peak_call(dbp)
    print "singlet", os.path.basename(sfname)
    sbp = get_tid2basep_from_ribomap_base(sfname, cds_range)
    speaks = batch_peak_call(sbp)
    dist_list = generate_dist_list(dpeaks, speaks, dstart, dend)
    return [ [tid, dpos, int(dbp[tid][dpos]), spos, int(sbp[tid][spos])] for tid,dpos,spos in dist_list]
    
def get_doublet_singlet_peak_codon_dist(dfname, sfname, cds_range, tseq, dstart=10, dend=40, seq_window=20):
    """
    return list: [transcript_id, doublet_loc, doublet_cnt, singlet_loc, singlet_cnt, codon ]
    default globals: dlen_min, dlen_max, peak_width, peak_min_cutoff
    """
    dcp = generate_codon_profile_from_rlen_hist(dfname, cds_range)
    dpeaks = batch_peak_call(dcp)
    print "singlet", os.path.basename(sfname)
    scp = get_tid2codonp_from_ribomap_base(sfname, cds_range)
    speaks = batch_peak_call(scp)
    dist_list = generate_dist_list(dpeaks, speaks, dstart, dend)
    return [ [tid, dpos, int(dcp[tid][dpos]), spos, int(scp[tid][spos]), tseq[tid][(spos-seq_window)*3:(spos+1+seq_window)*3]] for tid,dpos,spos in dist_list ]

def write_list_to_file(flist, header, fname):
    tf = open(fname,'w')
    tf.write(header)
    for f in flist:
        tf.write('\t'.join(map(str,f))+'\n')
    tf.close()

def main():
    """
    get a list of coupled peaks in singlet and doublet
    """
    cds_range = get_cds_range(cds_txt)    
    tseq = get_tseq(tfasta, cds_range)

    dstart = 10
    dend = 40
    odir = "../ds_cmp/"
    header = "transcript_id\tdoublet_loc\tdoublet_cnt\tsinglet_loc\tsinglet_cnt\tseq_in_20_codons\n"

    dist_list = get_doublet_singlet_peak_codon_dist(nchx_dfn, nchx_sfn, cds_range, tseq)
    fname = odir+"nchx.txt"
    write_list_to_file(dist_list, header, fname)
    
    #dist_list = get_doublet_singlet_peak_base_dist(chx_dfn, chx_sfn, cds_range)
    dist_list = get_doublet_singlet_peak_codon_dist(chx_dfn, chx_sfn, cds_range, tseq)
    fname = odir+"chx.txt"
    write_list_to_file(dist_list, header, fname)

    dist_list = get_doublet_singlet_peak_codon_dist(wt_dfn, wt_sfn, cds_range, tseq)
    fname = odir+"wt.txt"
    write_list_to_file(dist_list, header, fname)

    dist_list = get_doublet_singlet_peak_codon_dist(dom34_dfn, dom34_sfn, cds_range, tseq)
    fname = odir+"dom34.txt"
    write_list_to_file(dist_list, header, fname)

    # dlist = [ l[3]-l[1] for l in dist_list ]
    # if dstart > dend:
    #     bins = np.arange(dend-1, dstart+2)
    # else:
    #     bins = np.arange(dstart-1, dend+2)
    # ns,bins,patches = plt.hist(dlist,bins,alpha=0.5)
    # plt.savefig('test.pdf')

if __name__ == "__main__": main()



