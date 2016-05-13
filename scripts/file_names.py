#!/usr/bin/env python
cds_txt = "../ref/cds_range.txt"
tfasta = "../ref/protein_coding_100_filtered.fasta"
#==============================
# multi-mappings
#==============================
# Joel's data
nchx_dfn = "../rlen_hist/prelim/Lib-5-5-15_3_CACGAT_R1_nodup.hist"
nchx_sfn = "../ribomap/Lib-5-5-15_2_AGTTCC_R1_nodup.base"
nchx_nbc_dfn = "../rlen_hist/prelim/Lib-5-5-15_3_CACGAT_R1_nonempty.hist"
nchx_nbc_sfn = "../ribomap/Lib-5-5-15_2_AGTTCC_R1_nonempty.base"
chx_dfn = "../rlen_hist/prelim/Lib-5-5-15_7_ATCACG_R1_nodup.hist"
chx_sfn = "../ribomap/Lib-5-5-15_6_TCCCGA_R1_nodup.base"
chx_nbc_dfn = "../rlen_hist/prelim/Lib-5-5-15_7_ATCACG_R1_nonempty.hist"
chx_nbc_sfn = "../ribomap/Lib-5-5-15_6_TCCCGA_R1_nonempty.base"
# Guydosh's data
wt_dfn = "../rlen_hist/Guydosh14/Guydosh14_wt_disome.hist"
wt_sfn = "../ribomap/Guydosh14/Guydosh14_wt_CHX.base"
dom34_dfn = "../rlen_hist/Guydosh14/Guydosh14_dom34KO_disome.hist"
dom34_sfn = "../ribomap/Guydosh14/Guydosh14_dom34KO_CHX.base"
# significance test files
nchx_psfn = "../ds_cmp/nchx_singlet_pval.txt"
chx_psfn = "../ds_cmp/chx_singlet_pval.txt"
wt_psfn = "../ds_cmp/wt_singlet_pval.txt"
dom34_psfn = "../ds_cmp/dom34_singlet_pval.txt"
nchx_pdfn = "../ds_cmp/nchx_doublet_pval.txt"
chx_pdfn = "../ds_cmp/chx_doublet_pval.txt"
wt_pdfn = "../ds_cmp/wt_doublet_pval.txt"
dom34_pdfn = "../ds_cmp/dom34_doublet_pval.txt"

#==============================
# unique-mappings
#==============================
# Joel's data
nchx_dfn = "../uniquely_mapped/profiles/Lib-5-5-15_3_CACGAT_R1_nodup.hist"
nchx_sfn = "../uniquely_mapped/profiles/Lib-5-5-15_2_AGTTCC_R1_nodup_cds.base"
chx_dfn = "../uniquely_mapped/profiles/Lib-5-5-15_7_ATCACG_R1_nodup.hist"
chx_sfn = "../uniquely_mapped/profiles/Lib-5-5-15_6_TCCCGA_R1_nodup_cds.base"
# Guydosh's data
wt_dfn = "../uniquely_mapped/profiles/Guydosh14_wt_disome.hist"
wt_sfn = "../uniquely_mapped/profiles/Guydosh14_wt_CHX_cds.base"
dom34_dfn = "../uniquely_mapped/profiles/Guydosh14_dom34KO_disome.hist"
dom34_sfn = "../uniquely_mapped/profiles/Guydosh14_dom34KO_CHX_cds.base"
# significance test files
nchx_psfn = "../uniquely_mapped/jam_pvals/nchx_singlet_pval.txt"
chx_psfn = "../uniquely_mapped/jam_pvals/chx_singlet_pval.txt"
wt_psfn = "../uniquely_mapped/jam_pvals/wt_singlet_pval.txt"
dom34_psfn = "../uniquely_mapped/jam_pvals/dom34_singlet_pval.txt"
nchx_pdfn = "../uniquely_mapped/jam_pvals/nchx_doublet_pval.txt"
chx_pdfn = "../uniquely_mapped/jam_pvals/chx_doublet_pval.txt"
wt_pdfn = "../uniquely_mapped/jam_pvals/wt_doublet_pval.txt"
dom34_pdfn = "../uniquely_mapped/jam_pvals/dom34_doublet_pval.txt"

#==============================
# global params
#==============================
# doublet filter settings
dlen_min = 57
dlen_max = 62
# peak calling settings
peak_width=3
peak_min_cutoff=3
# significance test for collision
sd_distance = 13
window_size = 5
ds_distance = -(sd_distance + window_size)
sample_size = 1000000
