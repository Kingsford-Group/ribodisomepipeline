#!/usr/bin/env python
import os
def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)

def get_file_core(fname):
    istart = fname.rfind("/")
    iend = fname.rfind(".")
    return fname[istart:iend]
