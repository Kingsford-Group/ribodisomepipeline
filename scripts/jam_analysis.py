import scipy.spatial
import scipy.stats
import numpy as np
from compare_significance import *
from codon_peak import *

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['font.size'] = 15
rcParams['xtick.major.size'] = 5
rcParams['ytick.major.size'] = 5

def library_analysis_pipeline(psfname, pdfname, tseq, distance, window_size, peak_type, fn_short, oprfx):
    tid2ssig, tid2srest, tid2dsig, tid2drest, tid2jsig = prepare_peak_list(psfname, pdfname, distance, window_size, peak_type)
    codon_usage = get_codon_usage(tid2ssig.keys(), 1, 1, tseq, True)
    codon_peak = get_codon_peak_freq(tid2ssig, tseq)
    title="{0} singlet jam".format(fn_short)
    figprfx = "{0}_ssig_abs".format(oprfx)
    plot_codon_freq_bg_fg(codon_usage, codon_peak, title, figprfx)
    aa_usage = codon_cnt_to_aa_cnt(codon_usage)
    aa_peak = codon_cnt_to_aa_cnt(codon_peak)
    plot_aa_freq_bg_fg(aa_usage, aa_peak, title, figprfx)
    title="{0} singlet jam (normalized by aa)".format(fn_short)
    figprfx = "{0}_ssig_norm".format(oprfx)
    codon_nusage = normalize_cnt_by_aa(codon_usage)
    codon_npeak = normalize_cnt_by_aa(codon_peak)
    plot_codon_freq_bg_fg(codon_nusage, codon_npeak, title, figprfx)

    codon_usage = get_codon_usage(tid2jsig.keys(), 1, 1, tseq, True)
    codon_peak = get_codon_peak_freq(tid2jsig, tseq)
    title="{0} joint jam".format(fn_short)
    figprfx = "{0}_jsig_abs".format(oprfx)
    plot_codon_freq_bg_fg(codon_usage, codon_peak, title, figprfx)
    aa_usage = codon_cnt_to_aa_cnt(codon_usage)
    aa_peak = codon_cnt_to_aa_cnt(codon_peak)
    plot_aa_freq_bg_fg(aa_usage, aa_peak, title, figprfx)
    title="{0} joint jam (normalized by aa)".format(fn_short)
    figprfx = "{0}_jsig_norm".format(oprfx)
    codon_nusage = normalize_cnt_by_aa(codon_usage)
    codon_npeak = normalize_cnt_by_aa(codon_peak)
    plot_codon_freq_bg_fg(codon_nusage, codon_npeak, title, figprfx)

    codon_usage = get_codon_usage(tid2srest.keys(), 1, 1, tseq, True)
    codon_peak = get_codon_peak_freq(tid2srest, tseq)
    title="{0} singlet nonjam".format(fn_short)
    figprfx = "{0}_srest_abs".format(oprfx)
    plot_codon_freq_bg_fg(codon_usage, codon_peak, title, figprfx)
    aa_usage = codon_cnt_to_aa_cnt(codon_usage)
    aa_peak = codon_cnt_to_aa_cnt(codon_peak)
    plot_aa_freq_bg_fg(aa_usage, aa_peak, title, figprfx)
    title="{0} singlet nonjam (normalized by aa)".format(fn_short)
    figprfx = "{0}_srest_norm".format(oprfx)
    codon_nusage = normalize_cnt_by_aa(codon_usage)
    codon_npeak = normalize_cnt_by_aa(codon_peak)
    plot_codon_freq_bg_fg(codon_nusage, codon_npeak, title, figprfx)

def codon_dic_to_list(codon_list, codon_dic):
    return np.array([ codon_dic[codon] if codon in codon_dic else 0 for codon in codon_list ],dtype=float)

def get_tot_cnt(tid2peaks):
    return sum(map(len, tid2peaks.values()))

def get_expected_ratio(fg, bg):
    cnt_fg = get_tot_cnt(fg)
    cnt_bg = get_tot_cnt(bg)
    return float(cnt_fg)/(cnt_fg + cnt_bg)

def compute_chisquare(cnt_obs, cnt_exp):
    select = (cnt_exp != 0)
    vobs = cnt_obs[select]
    vexp = cnt_exp[select]
    chisq, p = scipy.stats.chisquare(vobs, vexp)
    return chisq, p, sum(select)-1

def rescale_vec(vadj, vref):
    return vadj*np.sum(vref)/np.sum(vadj)

def chisquare_codon_freq_from_pfname(pfname, tseq, figname):
    tid2sig, tid2rest = get_sig_peaks_from_file(pfname)
    codon_list = generate_cc_list(False)
    codon_peak = get_codon_peak_freq(tid2sig, tseq, False)
    codon_usage = get_codon_usage(tid2sig.keys(), 1, 1, tseq, False)
    sics_fg = codon_dic_to_list(codon_list, codon_peak)
    sics_bg = codon_dic_to_list(codon_list, codon_usage)
    codon_peak = get_codon_peak_freq(tid2rest, tseq, False)
    codon_usage = get_codon_usage(tid2rest.keys(), 1, 1, tseq, False)
    nonsics_fg = codon_dic_to_list(codon_list, codon_peak)
    nonsics_bg = codon_dic_to_list(codon_list, codon_usage)

    # rescale vectors since chi-square is sensitive to sample size
    sics_bg = rescale_vec(sics_bg, sics_fg)
    nonsics_fg = rescale_vec(nonsics_fg, sics_fg)
    nonsics_bg = rescale_vec(nonsics_bg, sics_fg)
    print sics_fg

    chisq, p, df = compute_chisquare(sics_fg, sics_bg)
    print "chisquare sics and bg {0:.0f} {1} {2}".format(chisq, p, df)
    chisq, p, df = compute_chisquare(nonsics_fg, nonsics_bg)
    print "chisquare nonsics and bg {0:.0f} {1} {2}".format(chisq, p, df)
    chisq, p, df = compute_chisquare(sics_fg, nonsics_fg)
    print "chisquare sics and nonsics {0:.0f} {1} {2}".format(chisq, p, df)
    chisq, p, df = compute_chisquare(nonsics_fg, sics_fg)
    print "chisquare nonsics and sics {0:.0f} {1} {2}".format(chisq, p, df)
    chisq, p, df = compute_chisquare(sics_bg, nonsics_bg)
    print "chisquare bg and bg {0:.0f} {1} {2}".format(chisq, p, df)

def codon_usage_sics_nonsics(pfname, tseq, figname):
    tid2sig, tid2rest = get_sig_peaks_from_file(pfname)
    codon_list = generate_cc_list()
    codon_usage = get_codon_usage(tid2rest.keys(), 1, 1, tseq, True)
    codon_peak = get_codon_peak_freq(tid2rest, tseq)
    nonsics_bg = codon_dic_to_list(codon_list, codon_usage)
    nonsics_fg = codon_dic_to_list(codon_list, codon_peak)
    print "nonsics to bg", scipy.spatial.distance.euclidean(nonsics_bg, nonsics_fg)
    codon_usage = get_codon_usage(tid2sig.keys(), 1, 1, tseq, True)
    codon_peak = get_codon_peak_freq(tid2sig, tseq)
    sics_bg = codon_dic_to_list(codon_list, codon_usage)
    sics_fg = codon_dic_to_list(codon_list, codon_peak)
    print "sics to bg", scipy.spatial.distance.euclidean(sics_bg, sics_fg)
    print "sics to nonsics", scipy.spatial.distance.euclidean(sics_fg, nonsics_fg)
    sics_ratio = np.ones(len(sics_bg))
    sics_ratio[sics_bg!=0] = sics_fg[sics_bg!=0]/sics_bg[sics_bg!=0]
    sics_ratio = np.log2(sics_ratio)
    print "plotting"
    aa2color = get_aa_colormap()
    c = [ aa2color[codon2aa[codon]] for codon in codon_list ]
    x = range(len(codon_list))
    plt.figure(figsize=(12,6))
    plt.bar(x, sics_ratio, color=c, edgecolor='white', align='center')
    plt.hlines([1,-1], x[0]-1, x[-1]+1, 'r')
    plt.xticks(x, codon_list, rotation='vertical', fontsize=12)
    plt.xlim((x[0]-1, x[-1]+1))
    plt.ylabel('codon usage ratio log2(fg/bg)')
    plt.tight_layout()
    plt.savefig(figname, bbox_inches='tight')
    plt.close()

if __name__ == "__main__":
    cds_range = get_cds_range(cds_txt)
    tseq = get_tseq(tfasta, cds_range)
    print "nchx"
    codon_usage_sics_nonsics(nchx_psfn, tseq, figure_dir+"nchx_sics_codon_ratio.pdf")
    chisquare_codon_freq_from_pfname(nchx_psfn, tseq, figure_dir+"nchx_sics_codon_freq.pdf")
    print "chx"
    codon_usage_sics_nonsics(chx_psfn, tseq, figure_dir+"chx_sics_codon_ratio.pdf")
    chisquare_codon_freq_from_pfname(chx_psfn, tseq, figure_dir+"chx_sics_codon_freq.pdf")
    print "wt"
    codon_usage_sics_nonsics(wt_psfn, tseq, figure_dir+"wt_sics_codon_ratio.pdf")
    chisquare_codon_freq_from_pfname(wt_psfn, tseq, figure_dir+"wt_sics_codon_freq.pdf")
    print "dom34"
    codon_usage_sics_nonsics(dom34_psfn, tseq, figure_dir+"dom34_sics_codon_ratio.pdf")
    chisquare_codon_freq_from_pfname(dom34_psfn, tseq, figure_dir+"dom34_sics_codon_freq.pdf")
    """
    distance = ds_distance
    peak_type = 'singlet'
    library_analysis_pipeline(nchx_psfn, nchx_pdfn, tseq, distance, window_size, peak_type, 'nchx', figure_dir+'nchx')
    library_analysis_pipeline(chx_psfn, chx_pdfn, tseq, distance, window_size, peak_type, 'chx', figure_dir+'chx')
    library_analysis_pipeline(wt_psfn, wt_pdfn, tseq, distance, window_size, peak_type, 'wt', figure_dir+'wt')
    library_analysis_pipeline(dom34_psfn, dom34_pdfn, tseq, distance, window_size, peak_type, 'dom34', figure_dir+'dom34')
    """
