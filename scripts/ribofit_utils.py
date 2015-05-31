#!/usr/bin/env python
import numpy as np
import scipy.stats
import exGaussian
from peak_cluster import threshold, identify_peaks, cluster_peaks

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['font.size'] = 20
rcParams['xtick.major.size'] = 5
rcParams['ytick.major.size'] = 5
cmap = matplotlib.cm.gist_rainbow

#=============================
# filter profiles
#=============================
def validate_Tuller(vec, threshold=1):
    # median > 1
    return np.median(vec) >= threshold

def validate_Pop(p):
    # <50% of positions on the CDS with at least 1 mapped
    # p is the profile values in the profile dic
    # p has keys: rid, rprofile, nprofile and mprofile
    return not ( np.mean(np.array(p['mprofile']) >= 1) < 0.5 or np.all(p['rprofile'] == 0) )

def validate_profile(vec):
    # > 50% cnts > 1
    return np.mean(vec>1)>0.5

#=============================
# profile transformation
#=============================
def merge_frames(p):
    """ merge off-frames to the closest frame0"""
    p_sum = np.array(p[::3])
    pplus = np.array(p[1::3])
    p_sum += pplus
    pplus = np.array([0]+p[2::3][:-1])
    p_sum += pplus
    return p_sum

def sum_frames(p):
    """ sum frame 1 and 2 to frame0 """
    return [ sum(p[i:i+3]) for i in xrange(0,len(p),3) ]

def in_frames(p):
    return p[::3]

def codonp_from_basep(p, cds_range, merge_func):
    print "converting base counts to codon counts..."
    pc = {}
    for rid in p:
        tid = p[rid]['tid']
        start, end = cds_range[tid]
        if (end-start)%3 != 0: 
            print "profile length not a duplicate of 3!", tid
            exit(1)
        pmr = merge_func(p[rid]['rprofile'][start:end])
        pmn = merge_func(p[rid]['nprofile'][start:end])
        pmm = p[rid]['mprofile'][start:end]
        pc[rid] = { 'tid' : tid,
                    'rprofile' : pmr,
                    'nprofile' : pmn,
                    'mprofile' : pmm}
    return pc

def merge_reps(p1,p2):
    p_sum = {}
    rid_list = list(set(p1.keys())|set(p2.keys()))
    for rid in rid_list:
        if rid not in p1:
            p_sum[rid] = p2[rid]
        else:
            p = np.array(p1[rid]['rprofile'])
            if rid in p2:
                p += np.array(p2[rid]['rprofile'])
            p_sum[rid] = { 'rprofile': p,
                           'tid': p1[rid]['tid'] }
    return p_sum

#=============================
# group counts by codon
#=============================
def include_loc(p, seq, loc, dic):
    codon = seq[loc*3 : (loc+1)*3]
    cnt = p[loc]
    append_nrc_to_codon(dic, codon, cnt)
    
def append_nrc_to_codon(dic, codon, nrc):
    dic.setdefault(codon, []).append(nrc)

def get_mean_rc(vec, threshold=1):
    return np.mean(vec[vec>threshold])

def group_nrc_by_codon(rc_list, p, tseq, rc_threshold=1, start=20, stop=-20):
    j = 0
    for rid in p:
        tid = p[rid]['tid']
        rprof = np.array(p[rid]['rprofile'])
        assert len(rprof)*3 == len(tseq[tid])
        # excluding first and last 20 codons of a transcript
        plen = len(rprof[start:stop])
        if plen==0: continue
        # discard transcript if median rc <= 1
        if not validate_Tuller(rprof[start:stop], rc_threshold): continue
        # compute average rc (by excluding rc <= 1)
        nrc = get_mean_rc(rprof[start:stop],rc_threshold)
        # group normalized rc by codon type
        for i in xrange(start, start+plen):
            if rprof[i] <= rc_threshold: continue
            codon = tseq[tid][i*3: (i+1)*3]
            cnt = rprof[i]/nrc
            append_nrc_to_codon(rc_list, codon, cnt)
        j += 1
    print "number of transcripts included: {0}".format(j)
    return rc_list


def group_np_tec_by_codon(p, tseq):
    j = 0
    jam_codon = {}
    jam_free_codon = {}
    for rid in p:
        tid = p[rid]['tid']
        rprof = np.array(p[rid]['rprofile'])
        assert len(rprof)*3 == len(tseq[tid])
        if not validate_profile(rprof): continue
        # only identify peaks
        peak = identify_peaks(rprof, threshold)
        prof = np.array(p[rid]['nprofile'])
        # group cnts by codons
        # skip first and last codon (start and stop)
        for c in xrange(1, len(prof)-1):
            if peak[c] == True :
                include_loc(prof, tseq[tid], c, jam_codon)
            elif rprof[c]>1 :
            #else:
                include_loc(prof, tseq[tid], c, jam_free_codon)
        j += 1
    print "total transcripts included: {0}".format(j)
    return jam_free_codon, jam_codon

def group_cnts_by_codon(jam_codon, jam_free_codon, p, tseq, rc_threshold=1, start=20, stop=-20, filter_peak=True, normalization="mRNA", validate_trans_func=validate_profile):
    j = 0
    for rid in p:
        tid = p[rid]['tid']
        rprof = np.array(p[rid]['rprofile'])
        assert len(rprof)*3 == len(tseq[tid])
        # to be consistent with Tuller, skip first and last 20 codons
        if stop == 0: stop = len(rprof)
        plen = len(rprof[start:stop])
        if plen==0: continue
        if not validate_trans_func(rprof[start:stop]): continue
        # identify peaks
        if filter_peak == True:
            peak = identify_peaks(rprof, threshold)
        if normalization == "mRNA":
            prof = np.array(p[rid]['nprofile'])
        elif normalization == "nrc":
            nrc = get_mean_rc(rprof[start:stop], rc_threshold)
            prof = rprof/nrc
        else:
            print "unrecoginized normalization method {0}!".format(normalization)
            exit(1)
        # group cnts by codons
        for c in xrange(start, start+plen):
            if filter_peak == True:
                if peak[c] == True:
                    include_loc(prof, tseq[tid], c, jam_codon)
                elif rprof[c] > rc_threshold :
                    include_loc(prof, tseq[tid], c, jam_free_codon)
            elif rprof[c]> rc_threshold :
                include_loc(prof, tseq[tid], c, jam_free_codon)
        j += 1
    print "total transcripts included: {0}".format(j)
    return jam_free_codon, jam_codon

def group_npregion_tec_by_codon(p, tseq, ribo_len):
    j = 0
    jam_codon = {}
    jam_free_codon = {}
    for rid in p:
        tid = p[rid]['tid']
        rprof = np.array(p[rid]['rprofile'])
        assert len(rprof)*3 == len(tseq[tid])
        if not validate_profile(rprof): continue
        # identify peak regions
        peak = identify_peaks(rprof, threshold)
        seg_list = cluster_peaks(peak,ribo_len*2,ribo_len)
        # group cnts by codons
        prof = np.array(p[rid]['nprofile'])
        # list empty: no peak-rich region
        if not seg_list:
            for c in xrange(len(prof)):
                include_loc(prof, tseq[tid], c, jam_free_codon)
        # list not empty, separate two classes of regions
        else:
            # 1st jam-free seg before 1st jam seg
            for c in xrange(seg_list[0][0]):
                include_loc(prof, tseq[tid], c, jam_free_codon)
            # go through each jam seg and jam-free seg after it
            for i in xrange(len(seg_list)-1):
                for c in xrange(seg_list[i][0], seg_list[i][1]+1):
                    include_loc(prof, tseq[tid], c, jam_codon)
                for c in xrange(seg_list[i][1]+1, seg_list[i+1][0]):
                    include_loc(prof, tseq[tid], c, jam_free_codon)
            # last jam seg and jam-free seg after
            for c in xrange(seg_list[-1][0], seg_list[-1][1]+1):
                include_loc(prof, tseq[tid], c, jam_codon)
            for c in xrange(seg_list[-1][1]+1, len(prof)):
                include_loc(prof, tseq[tid], c, jam_free_codon)
        j += 1
    print "total transcripts included: {0}".format(j)
    return jam_free_codon, jam_codon
    
#=============================
# model fitting
#=============================
def fit_exGaussian(rc_list):
    params = {}
    for codon in rc_list:
        nrc_list = np.array(rc_list[codon])
        # percentile cutoff to garantee valid analytical solution
        threshold = np.percentile(nrc_list, 95)
        #threshold = 3
        mu, sigma, lamb = exGaussian.fit(nrc_list[nrc_list<threshold])
        params[codon] = (mu, sigma, lamb)
    return params
    
def fit_lognormal(rc_list):
    params = {}
    for codon in rc_list:
        shape, loc, scale = scipy.stats.lognorm.fit(rc_list[codon])
        params[codon] = (shape, loc, scale)
    return params

def fit_expon(rc_list):
    params = {}
    for codon in rc_list:
        loc, scale = scipy.stats.expon.fit(rc_list[codon])
        params[codon] = (loc, scale)
    return params

#=============================
# evaluate fitting
#=============================
def exclude_keys(dic, key_list):
    return { k:v for k,v in dic.iteritems() if k not in key_list }

def negloglikelihood(params, data, pdf_func):
    pdf = pdf_func(data, params)
    logpdf = np.log(pdf[pdf!=0])
    nllk = -np.sum(logpdf)
    dcnt = len(logpdf)
    return nllk, dcnt

def AIC(params, data, pdf_func):
    k = len(params)
    nllk, dcnt = negloglikelihood(params, data, pdf_func)
    aic = 2*k + 2*nllk
    return aic

def tot_aic(params_list, rc_list, pdf_func):
    aic = 0
    for codon in params_list:
        if codon not in rc_list: continue
        aic_cur = AIC(params_list[codon], rc_list[codon], pdf_func)
        aic += aic_cur
    return aic

def BIC(params, data, pdf_func):
    nllk = 0
    nparam = 0
    ndata = 0
    for codon in params:
        k = len(params[codon])
        nparam += k
        if codon not in data: continue
        nllk_cur, dcnt = negloglikelihood(params[codon], data[codon], pdf_func)
        nllk += nllk_cur
        ndata += dcnt
    return 2 * nllk + nparam * np.log(ndata), nllk, nparam, ndata

def sum_ksd(params, data, cdf_func):
    ks_list = []
    for codon in params:
        if codon not in data: continue
        data_c = data[codon]
        param_c = params[codon]
        d, p = scipy.stats.kstest(data_c, cdf_func, param_c)
        ks_list.append(d)
    return np.sum(ks_list)

def evaluate_fit_pipeline(jam_free_codon, jam_codon, tAI, xmin, xmax, fn_prefix):
    clist = tAI.keys()
    tai_list = np.array([tAI[c] for c in clist])

    print "fit footprint to exponentially modified Gaussian distributions"
    params = fit_exGaussian(jam_free_codon)
    plot_codon_hist_fitting(jam_free_codon, params, exGaussian.pdf, jam_codon, xmin, xmax, fn_prefix+"_emg")
    bic, nllk, nparam, ndata = BIC(params, jam_free_codon, exGaussian.pdf)
    print "BIC: {0:.2f} nllk: {1:.2f} nparam: {2} ndata {3}".format(bic, nllk, nparam, ndata)
    if len(jam_codon)!=0:
        bicj, nllkj, nparamj, ndataj = BIC(params, jam_codon, exGaussian.pdf)
        print "with peaks: BIC: {0:.2f} nllk: {1:.2f} nparam: {2} ndata {3}".format(2*(nllk+nllkj)+nparam*(ndata+ndataj), nllk+nllkj, nparam, ndata+ndataj)
    aic = tot_aic(params, jam_free_codon, exGaussian.pdf)
    print "AIC: {0:.2f} data points {1} average AIC: {2:.2f}".format(aic, ndata, aic/ndata)
    print "sum KS: {0:.2f}".format(sum_ksd(params, jam_free_codon, exGaussian.cdf))
    print "compare parameters with tAI"
    mu = [ params[c][0] for c in clist ]
    scatter_tai(tai_list, mu, "mu", fn_prefix+"_emg_mu")
    sigma = [ params[c][1] for c in clist ]
    scatter_tai(tai_list, sigma, "sigma", fn_prefix+"_emg_sigma")
    lamb = [ params[c][2] for c in clist ]
    scatter_tai(tai_list, lamb, "1/lambda", fn_prefix+"_emg_inverse_lambda")
    
    print "fit footprint to lognormal distributions"
    params = fit_lognormal(jam_free_codon)
    plot_codon_hist_fitting(jam_free_codon, params, lognormal_pdf, jam_codon, xmin, xmax, fn_prefix+"_lognorm")
    bic, nllk, nparam, ndata = BIC(params, jam_free_codon, lognormal_pdf)
    print "BIC: {0:.2f} nllk: {1:.2f} nparam: {2} ndata {3}".format(bic, nllk, nparam, ndata)
    if len(jam_codon)!=0:
        bicj, nllkj, nparamj, ndataj = BIC(params, jam_codon, lognormal_pdf)
        print "with peaks: BIC: {0:.2f} nllk: {1:.2f} nparam: {2} ndata {3}".format(2*(nllk+nllkj)+nparam*(ndata+ndataj), nllk+nllkj, nparam, ndata+ndataj)
    aic = tot_aic(params, jam_free_codon, lognormal_pdf)
    print "AIC: {0:.2f} data points {1} average AIC: {2:.2f}".format(aic, ndata, aic/ndata)
    print "sum KS: {0:.2f}".format(sum_ksd(params, jam_free_codon, scipy.stats.lognorm.cdf))
    print "compare parameters with tAI"
    sigma = [ params[c][0] for c in clist ]
    scatter_tai(tai_list, sigma, "sigma", fn_prefix+"_lognorm_sigma")
    mu = [ np.log(params[c][2]) for c in clist ]
    scatter_tai(tai_list, mu, "mu", fn_prefix+"_lognorm_mu")
    median = np.exp(np.array(mu))
    scatter_tai(tai_list, median, "median", fn_prefix+"_lognorm_median")
    var = np.square(np.array(sigma))
    skew = (np.exp(var)+2)*(np.sqrt(np.exp(var)-1))
    scatter_tai(tai_list, skew, "skewness", fn_prefix+"_lognorm_skew")

    print "fit footprint to exponential distributions"    
    params = fit_expon(jam_free_codon)
    plot_codon_hist_fitting(jam_free_codon, params, expon_pdf, jam_codon, xmin, xmax, fn_prefix+"_exp")
    bic, nllk, nparam, ndata = BIC(params, jam_free_codon, expon_pdf)
    print "BIC: {0:.2f} nllk: {1:.2f} nparam: {2} ndata {3}".format(bic, nllk, nparam, ndata)
    if len(jam_codon)!=0:
        bicj, nllkj, nparamj, ndataj = BIC(params, jam_codon, expon_pdf)
        print "with peaks: BIC: {0:.2f} nllk: {1:.2f} nparam: {2} ndata {3}".format(2*(nllk+nllkj)+nparam*(ndata+ndataj), nllk+nllkj, nparam, ndata+ndataj)
    aic = tot_aic(params, jam_free_codon, expon_pdf)
    print "AIC: {0:.2f} data points {1} average AIC: {2:.2f}".format(aic, ndata, aic/ndata)
    print "sum KS: {0:.2f}".format(sum_ksd(params, jam_free_codon, scipy.stats.expon.cdf))    
    lamb = [ params[c][1] for c in clist ]
    scatter_tai(tai_list, lamb, "1/lambda", fn_prefix+"_exp_inverse_lambda")

#=============================
# plot
#=============================
def lognormal_pdf(x, params):
    """ wrapper function for plotting"""
    shape, loc, scale = params
    return scipy.stats.lognorm.pdf(x,shape,loc, scale)

def expon_pdf(x, params):
    loc, scale = params
    return scipy.stats.expon.pdf(x,loc, scale)

def plot_codon_hist_fitting(rc_list, params, pdf_func, rc_jam_list=None, xmin=1e-7, xmax=3, fn_prefix="fp"):
    print "plot histograms"
    fig = plt.figure(facecolor='white', figsize=(16, 16))
    i = 0
    clist = sorted(rc_list.keys())
    for codon in clist:
        nrc_list = np.array(rc_list[codon])
        i += 1
        ax = fig.add_subplot(8, 8, i)
        ax.set_frame_on(False)
        ax.get_xaxis().tick_bottom()
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().tick_left()
        ax.get_yaxis().set_visible(False)
        ns,bins,patches = plt.hist(nrc_list, np.linspace(xmin,xmax,50), normed = True, histtype='stepfilled', log=False, color = 'c', edgecolor='c', alpha=0.5)
        ymax = max(ns)
        if rc_jam_list!=None and codon in rc_jam_list and len(rc_jam_list[codon])>10:
            ns,bins,patches = plt.hist(rc_jam_list[codon], np.linspace(xmin,xmax,50), normed = True, histtype='stepfilled', log=False, color='b', edgecolor='b', hatch='/', alpha=0.3)
            ymax = max(ns+[ymax])
        x = np.linspace(xmin, xmax, 10000)
        y = pdf_func(x, params[codon])
        plt.plot(x,y,'r-',lw=2, alpha=0.3)
        ax.set_title(codon,fontsize=10)
        ax.set_xlim((0,xmax))
        ax.set_ylim((0,ymax+0.1))
    plt.savefig(fn_prefix+"_codon_hist.png", bbox_inches='tight')
    plt.close()
    #plt.show()

def scatter_tai(tai_list, vals, val_name, fn_prefix):
    plt.figure()
    plt.plot(tai_list, vals, 'bo', alpha=0.5)
    plt.xlabel('tAI')
    plt.ylabel("{0}".format(val_name))
    pr = scipy.stats.pearsonr(tai_list,vals)[0]
    sr = scipy.stats.spearmanr(tai_list, vals)[0]
    linear_fit = np.polyfit(tai_list, vals,1)
    fit_fn = np.poly1d(linear_fit)
    plt.plot(tai_list, fit_fn(tai_list),c='r',linewidth=2, alpha=0.5)
    plt.title("pearson r = {0:.2f} spearman r = {1:.2f}".format(pr, sr))
    plt.savefig(fn_prefix+"_scatter.pdf", bbox_inches='tight')
    #plt.show()
    print "{0}: pearsonr: {1:.2f} spearmanr {2:.2f}".format(val_name, pr, sr)
    
#=============================
# IO
#=============================
def write_params(fn, params):
    tf = open(fn,'w')
    for codon in sorted(params.keys()):
        txt = [codon]+[str(i) for i in params[codon]]
        tf.write('\t'.join(txt)+'\n')
    tf.close()

def read_params(fn, line_sep='\n'):
    params = {}
    tf = open(fn)
    txt = tf.read()
    lines = txt.split(line_sep)
    for line in lines:
        words = line.split()
        params[words[0]] = map(float, words[1:])
    tf.close()
    return params
