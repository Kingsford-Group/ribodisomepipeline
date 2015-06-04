#!/usr/bin/env python
import sys
import numpy as np
import scipy.stats
from ribomap_result_parser import *
from ribofit_utils import *
from io_utils import *

def main():
    if len(sys.argv)!=5:
        print "Usage: python codon_count_hist.py cds_range.txt ref.fasta ribomap.base odir"
        exit(1)

    cds_txt = sys.argv[1]
    tfasta = sys.argv[2]
    pfn = sys.argv[3]
    odir = sys.argv[4]
    cds_range = get_cds_range(cds_txt)
    tseq = get_tseq(tfasta, cds_range)
    pb=parse_estimated_profile(pfn)
    rc_threshold = 1
    start=20
    stop=-20
    
    # tAI = get_rates("absolute_adaptiveness_SCer.txt")
    # stop_codons = ["TAA", "TGA", "TAG"]
    # tAI = exclude_keys(tAI, stop_codons)

    ensure_dir(odir)
    print "sum counts from all three frames"
    pc = codonp_from_basep(pb, cds_range, sum_frames)

    print "Tuller's filtering and normalization scheme"
    jam_codon = {}
    jam_free_codon = {}
    jam_free_codon, jam_codon = group_cnts_by_codon(jam_codon, jam_free_codon, pc, tseq,
                                                    rc_threshold=rc_threshold, start=start, stop=stop,
                                                    filter_peak=False, normalization="nrc",
                                                    validate_trans_func=validate_profile)

    print "fitting lognormal"
    params = fit_lognormal(jam_free_codon)
    fn_prefix = odir+get_file_core(pfn)+"_lognormal"
    plot_codon_hist_fitting(jam_free_codon, params, lognormal_pdf, fn_prefix=fn_prefix)

    print "fitting exponential"
    jam_codon = {}
    jam_free_codon = {}
    jam_free_codon, jam_codon = group_cnts_by_codon(jam_codon, jam_free_codon, pc, tseq,
                                                    rc_threshold=-1, start=start, stop=stop,
                                                    filter_peak=True, normalization="mRNA",
                                                    validate_trans_func=validate_profile)
    params = fit_expon(jam_free_codon)
    fn_prefix = odir+get_file_core(pfn)+"_exp"
    plot_codon_hist_fitting(jam_free_codon, params, expon_pdf, jam_codon, fn_prefix=fn_prefix)

if __name__ == "__main__": main()
    


