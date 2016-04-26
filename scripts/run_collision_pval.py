#!/usr/bin/env python
import os
from file_names import *

ofname = "../ds_cmp/nchx_singlet_pval.txt"
cmd = "python significant_doublet_count.py {0} {1} {2} {3} {4} {5}&".format(nchx_dfn, nchx_sfn, ofname, sd_distance, window_size, 'singlet')
os.system(cmd)

ofname = "../ds_cmp/nchx_doublet_pval.txt"
cmd = "python significant_doublet_count.py {0} {1} {2} {3} {4} {5}&".format(nchx_dfn, nchx_sfn, ofname, ds_distance, window_size, 'doublet')
os.system(cmd)

ofname = "../ds_cmp/chx_singlet_pval.txt"
cmd = "python significant_doublet_count.py {0} {1} {2} {3} {4} {5}&".format(chx_dfn, chx_sfn, ofname, sd_distance, window_size, 'singlet')
os.system(cmd)

ofname = "../ds_cmp/chx_doublet_pval.txt"
cmd = "python significant_doublet_count.py {0} {1} {2} {3} {4} {5}&".format(chx_dfn, chx_sfn, ofname, ds_distance, window_size, 'doublet')
os.system(cmd)

ofname = "../ds_cmp/wt_singlet_pval.txt"
cmd = "python significant_doublet_count.py {0} {1} {2} {3} {4} {5}&".format(wt_dfn, wt_sfn, ofname, sd_distance, window_size, 'singlet')
os.system(cmd)

ofname = "../ds_cmp/wt_doublet_pval.txt"
cmd = "python significant_doublet_count.py {0} {1} {2} {3} {4} {5}&".format(wt_dfn, wt_sfn, ofname, ds_distance, window_size, 'doublet')
os.system(cmd)

ofname = "../ds_cmp/dom34_singlet_pval.txt"
cmd = "python significant_doublet_count.py {0} {1} {2} {3} {4} {5}&".format(dom34_dfn, dom34_sfn, ofname, sd_distance, window_size, 'singlet')
os.system(cmd)

ofname = "../ds_cmp/dom34_doublet_pval.txt"
cmd = "python significant_doublet_count.py {0} {1} {2} {3} {4} {5}&".format(dom34_dfn, dom34_sfn, ofname, ds_distance, window_size, 'doublet')
os.system(cmd)
