from compare_significance import *
from codon_peak import *

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

if __name__ == "__main__":
    cds_range = get_cds_range(cds_txt)
    tseq = get_tseq(tfasta, cds_range)
    distance = ds_distance
    peak_type = 'singlet'
    core = ['nchx', 'chx', 'wt', 'dom34']
    psfname = [ "../ds_cmp/{0}_singlet_pval.txt".format(c) for c in core ]
    pdfname = [ "../ds_cmp/{0}_doublet_pval.txt".format(c) for c in core ]
    oprfx = "../figures/codon_peak_cnt/"
    for i in range(len(core)):
        library_analysis_pipeline(psfname[i], pdfname[i], tseq, distance, window_size, peak_type, core[i], oprfx+core[i])
