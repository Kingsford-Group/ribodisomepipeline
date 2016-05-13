#!/usr/bin/env python
import os
import sys
import numpy as np
from Bio import SeqIO

#==============================
# misc
#==============================
def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)

def get_file_core(fname):
    istart = fname.rfind("/")
    iend = fname.rfind(".")
    return fname[istart:iend]

#==============================
# rlen.hist parser
#==============================
def parse_rlen_hist(fname):
    print "parsing footprint counts grouped by read length"
    tlist = {}
    transcript = {}
    tf = open(fname)
    line = tf.readline()
    while line:
        if line.startswith("refID: "):
            if transcript:
                rid = transcript["id"]
                sys.stdout.write("processed transcript {0}.\t\r".format(rid))
                sys.stdout.flush()
                tlist[rid] = transcript.copy()
                transcript.clear()
            rid = int(line.lstrip("refID: ").rstrip("\n"))
            transcript["id"] = rid
        elif line.startswith("tid: "):
            tid = line.lstrip("tid: ").rstrip("\n")
            transcript["tid"] = tid
            transcript["prof"] = {}
        # profile len grouped by read length
        # format: read len: pos count, ...
        # [pos count] is the adjacency list way to save space
        else:
            read_len, pc_str = line.rstrip("\n").split(": ")
            read_len = int(read_len)
            prof = []
            pc_pairs = pc_str.rstrip(", ").split(", ")
            for pc in pc_pairs:
                pos, cnt = pc.split(" ")
                prof.append([int(pos), float(cnt)])
            transcript["prof"][read_len] = prof
        line = tf.readline()
    if transcript:
        rid = transcript["id"]
        tlist[rid] = transcript
    tf.close()
    sys.stdout.write('\n')
    return tlist

#------------------------------
# stats on read length
# and frame distribution
#------------------------------
def print_stats(tlist, cds_range):
    print "get read frame stats"
    fcnt = [0]*3
    rlmin = 100
    rlmax = 0
    for rid in tlist:
        for rlen in tlist[rid]['prof']:
            if rlen < rlmin:
                rlmin = rlen
            if rlen > rlmax:
                rlmax = rlen
            tid = tlist[rid]['tid']
            start, stop = cds_range[tid]
            for pos, cnt in tlist[rid]['prof'][rlen]:
                fcnt[(pos-start)%3] += cnt
    tcnt = float(sum(fcnt))
    print "total footprints: {0:.0f}".format(tcnt)
    print "frame distribution: {0:.2%} {1:.2%} {2:.2%}".format(fcnt[0]/tcnt, fcnt[1]/tcnt, fcnt[2]/tcnt)
    print "read length range: {0}-{1}".format(rlmin, rlmax)
    return rlmin, rlmax

#==============================
# deblur output parser
#==============================
def parse_deblur_profile(fname):
    print "start parsing deblurred profile..."
    pcds = {}
    tf = open(fname)
    line = tf.readline()
    while line:
        if line.startswith("tid: "):
            tid = line.lstrip("tid: ").rstrip("\n")
            line = tf.readline()
            prof = map(float, line.rstrip().split())
            pcds[tid] = np.array(prof)
        line = tf.readline()
    tf.close()
    return pcds

#------------------------------
# frame distribution stats
#------------------------------
def print_frame_stats(pribo):
    f = np.zeros(3)
    for tid, prof in pribo.iteritems():
        for i in xrange(3):
            f[i] += np.sum(prof[i::3])
    f /= np.sum(f)
    print "{0:.2%} {1:.2%} {2:.2%}".format(f[0], f[1], f[2])

#==============================
# ribomap result parser
#==============================
def parse_estimated_profile(fname):
    print "start parsing ribomap profile..."
    pest = {}
    transcript = {}
    tf = open(fname)
    line = tf.readline()
    while line:
        if line.startswith("refID: "):
            if transcript:
                rid = transcript["id"]
                pest[rid] = transcript.copy()
                transcript.clear()
            rid = int(line.lstrip("refID: ").rstrip("\n"))
            sys.stdout.write("processed transcript {0}.\t\r".format(rid))
            sys.stdout.flush()
            transcript["id"] = rid
        elif line.startswith("tid: "):
            tid = line.lstrip("tid: ").rstrip("\n")
            transcript["tid"] = tid
        elif line.startswith("ribo profile: "):
            rprofile = map(float, line.lstrip("ribo profile: ").rstrip("\n").split())
            transcript["rprofile"] = rprofile
        elif line.startswith("normalized ribo profile: "):
            nprofile = map(float, line.lstrip("normalized ribo profile: ").rstrip("\n").split())
            transcript["nprofile"] = nprofile
        elif line.startswith("mRNA profile: "):
            mprofile = map(float, line.lstrip("mRNA profile: ").rstrip("\n").split())
            transcript["mprofile"] = mprofile
        else:
            pass
        line = tf.readline()
    if transcript:
        rid = transcript["id"]
        pest[rid] = transcript
    tf.close()
    sys.stdout.write('\n')
    return pest

#==============================
# cds range parser
#==============================
def get_cds_range(fn):
    print "getting cds range.."
    ifile = open(fn, 'r')
    cds_range = {}
    for line in ifile:
        tid, start, stop = line.strip().split()
        cds_range[tid] = (int(start), int(stop))
    ifile.close()
    return cds_range

#==============================
# fasta seq parser
#==============================
def get_tseq(fn, cds_range):
    print "getting transcript sequences..."
    ifile = open(fn, 'r')
    tseq = {}
    for rec in SeqIO.parse(ifile, "fasta"):
        tid = rec.id.split()[0]
        start, stop = cds_range[tid]
        seq = str(rec.seq[start:stop])
        tseq[tid] = seq
    ifile.close()
    return tseq

def get_pep_seq(fn):
    print "getting peptide sequence..."
    ifile = open(fn, 'r')
    pep_seq = {}
    for rec in SeqIO.parse(ifile, "fasta"):
        tid = rec.id.split()[0]
        pseq = str(rec.seq)
        pep_seq[tid] = pseq
    ifile.close()
    return pep_seq

