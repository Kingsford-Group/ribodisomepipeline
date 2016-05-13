#!/usr/bin/env python
import os
from file_names import *

multimap='False'

cmd = "python significant_doublet_count.py {0} {1} {2} {3} {4} {5} {6}&".format(nchx_dfn, nchx_sfn, nchx_psfn, sd_distance, window_size, 'singlet', multimap)
os.system(cmd)

cmd = "python significant_doublet_count.py {0} {1} {2} {3} {4} {5} {6}&".format(nchx_dfn, nchx_sfn, nchx_pdfn, ds_distance, window_size, 'doublet', multimap)
os.system(cmd)

cmd = "python significant_doublet_count.py {0} {1} {2} {3} {4} {5} {6}&".format(chx_dfn, chx_sfn, chx_psfn, sd_distance, window_size, 'singlet', multimap)
os.system(cmd)

cmd = "python significant_doublet_count.py {0} {1} {2} {3} {4} {5} {6}&".format(chx_dfn, chx_sfn, chx_pdfn, ds_distance, window_size, 'doublet', multimap)
os.system(cmd)

cmd = "python significant_doublet_count.py {0} {1} {2} {3} {4} {5} {6}&".format(wt_dfn, wt_sfn, wt_psfn, sd_distance, window_size, 'singlet', multimap)
os.system(cmd)

cmd = "python significant_doublet_count.py {0} {1} {2} {3} {4} {5} {6}&".format(wt_dfn, wt_sfn, wt_pdfn, ds_distance, window_size, 'doublet', multimap)
os.system(cmd)

cmd = "python significant_doublet_count.py {0} {1} {2} {3} {4} {5} {6}&".format(dom34_dfn, dom34_sfn, dom34_psfn, sd_distance, window_size, 'singlet', multimap)
os.system(cmd)

cmd = "python significant_doublet_count.py {0} {1} {2} {3} {4} {5} {6}&".format(dom34_dfn, dom34_sfn, dom34_pdfn, ds_distance, window_size, 'doublet', multimap)
os.system(cmd)
