#!/usr/bin/env python
from io_utils import get_tseq, get_cds_range
from file_names import *
from significant_doublet_count import read_pvals_from_file, get_window_cnt
from compare_significance import prepare_peak_list, max_upstream
from doublet_profile import generate_codon_profile_from_rlen_hist, get_tid2codonp_from_ribomap_base, get_codon_profile_from_deblur_profile

def get_region(seq, loc, offset, region_len):
    """
    extract sequence near a location with specific offset and length
    location unit: codon
    """
    istart = loc + offset
    if istart < 0: istart = 0
    iend = istart + region_len
    if iend > len(seq)/3: iend = len(seq)/3
    if istart == iend: return ''
    else: return seq[istart*3 : iend*3]

def get_batch_region(tid2loci, tseq, offset, region_len):
    tid2seq = {}
    for tid, loc_list in tid2loci.iteritems():
        tid2seq[tid] = {}
        for loc in loc_list:
            seq = get_region(tseq[tid], loc, offset, region_len)
            if seq != '':
                tid2seq[tid][loc] = seq
    return tid2seq

def get_seq_range(seq, istart, iend):
    if istart < 0: istart = 0
    if iend > len(seq)/3: iend = len(seq)/3
    if istart == iend: return ''
    else: return seq[istart*3 : iend*3]

def get_batch_merged_region(tid2loci, tseq, offset, region_len):
    tid2seq = {}
    for tid, loc_list in tid2loci.iteritems():
        tid2seq[tid] = {}
        loc_org = sorted(loc_list)
        loc_str = str(loc_org[0])
        istart = loc_org[0] + offset
        iend = istart + region_len
        for i in range(1,len(loc_org)):
            idelta = loc_org[i] - loc_org[i-1]
            if idelta < region_len:
                iend += idelta
                loc_str += "_" + str(loc_org[i])
                
            else:
                seq = get_seq_range(tseq[tid], istart, iend)
                if seq != '':
                    tid2seq[tid][loc_str] = seq
                loc_str = str(loc_org[i])
                istart = loc_org[i] + offset
                iend = istart + region_len
        if loc_str not in tid2seq[tid]:
            seq = get_seq_range(tseq[tid], istart, iend)
            if seq != '':
                tid2seq[tid][loc_str] = seq
    return tid2seq

def write_region(ofa, tid2seq):
    tf = open(ofa, 'w')
    for tid, seq_list in tid2seq.iteritems():
        for loc, seq in seq_list.iteritems():
            seq_id = ">{0}_{1}".format(tid, loc)
            tf.write("{0}\n{1}\n".format(seq_id,seq))
    tf.close()

def get_cnt_info(tid2loci, dcp, scp):
    tid2cinfo = {}
    for tid, loci_list in tid2loci.iteritems():
        tid2cinfo[tid] = {}
        tid2cinfo[tid]['len'] = len(scp[tid])
        tid2cinfo[tid]['sld'] = sum(scp[tid])
        tid2cinfo[tid]['dld'] = sum(dcp[tid])
        for loc in loci_list:
            scnt = scp[tid][loc]
            dmax = max_upstream(loc, dcp[tid])
            dsum = get_window_cnt(loc, dcp[tid])
            tid2cinfo[tid][loc] = [ scnt, dmax, dsum ]
    return tid2cinfo
            
def write_cnt_seq(ofname, tid2seqs, tid2cinfo, tid2jsig):
    header = "transcript_id\ttranscript_len\tsinglet_loads\tdoublet_loads\tsinglet_jam_loc\tsinglet_cnt\tmax_doublet_cnt\tsum_doublet_cnt\tjoint_jam\tseq\n"
    tf = open(ofname, 'w')
    tf.write(header)
    for tid in tid2seqs:
        for loc in tid2seqs[tid]:
            if tid not in tid2jsig:
                js = 'False'
            elif loc in tid2jsig[tid]:
                js = 'True'
            else:
                js = 'False'
            line = "{0}\t{1}\t{2:.0f}\t{3:.0f}\t{4:.0f}\t{5:.0f}\t{6:.0f}\t{7:.0f}\t{8}\t{9}\n".format(tid, tid2cinfo[tid]['len'], tid2cinfo[tid]['sld'], tid2cinfo[tid]['dld'], loc, tid2cinfo[tid][loc][0], tid2cinfo[tid][loc][1], tid2cinfo[tid][loc][2], js, tid2seqs[tid][loc])
            tf.write(line)
    tf.close()

def jam_info_pipeline(psfname, pdfname, sfname, dfname, distance, window_size, peak_type, offset, region_len, ofname, multimap):
    cds_range = get_cds_range(cds_txt)
    tseq = get_tseq(tfasta, cds_range)
    tid2ssig, tid2srest, tid2dsig, tid2drest, tid2jsig = prepare_peak_list(psfname, pdfname, distance, window_size, peak_type)
    tid2seqs = get_batch_region(tid2ssig, tseq, offset, region_len)
    dcp = generate_codon_profile_from_rlen_hist(dfname, cds_range)
    if multimap == True:
        scp = get_tid2codonp_from_ribomap_base(sfname, cds_range)
    else:
        scp = get_codon_profile_from_deblur_profile(sfname)
    tid2cinfo = get_cnt_info(tid2ssig, dcp, scp)
    write_cnt_seq(ofname, tid2seqs, tid2cinfo, tid2jsig)

def jam_fasta_pipeline(psfname, pdfname, distance, window_size, peak_type, offset, region_len, ofa):
    cds_range = get_cds_range(cds_txt)
    tseq = get_tseq(tfasta, cds_range)
    tid2ssig, tid2srest, tid2dsig, tid2drest, tid2jsig = prepare_peak_list(psfname, pdfname, distance, window_size, peak_type)
    tid2seqs = get_batch_merged_region(tid2ssig, tseq, offset, region_len)
    write_region(ofa, tid2seqs)

if __name__ == "__main__":
    distance = ds_distance
    peak_type = 'singlet'
    offset = -19
    region_len = 20

    jam_fasta_pipeline(nchx_psfn, nchx_pdfn, distance, window_size, peak_type, offset, region_len, nchx_fa)
    jam_fasta_pipeline(chx_psfn, chx_pdfn, distance, window_size, peak_type, offset, region_len, chx_fa)
    jam_fasta_pipeline(wt_psfn, wt_pdfn, distance, window_size, peak_type, offset, region_len, wt_fa)
    jam_fasta_pipeline(dom34_psfn, dom34_pdfn, distance, window_size, peak_type, offset, region_len, dom34_fa)

    jam_info_pipeline(nchx_psfn, nchx_pdfn, nchx_sfn, nchx_dfn, distance, window_size, peak_type, offset, region_len, nchx_info, multimap)
    jam_info_pipeline(chx_psfn, chx_pdfn, chx_sfn, chx_dfn, distance, window_size, peak_type, offset, region_len, chx_info, multimap)
    jam_info_pipeline(wt_psfn, wt_pdfn, wt_sfn, wt_dfn, distance, window_size, peak_type, offset, region_len, wt_info, multimap)
    jam_info_pipeline(dom34_psfn, dom34_pdfn, dom34_sfn, dom34_dfn, distance, window_size, peak_type, offset, region_len, dom34_info, multimap)

