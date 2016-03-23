#!/usr/bin/env python
from Bio import SeqIO
import sys
from collections import defaultdict

def parse_rlen_hist(fname):
    print "parsing footprint counts grouped by read length"
    tlist = {}
    transcript = {}
    tf = open(fname)
    line = tf.readline()
    while line:
        if line.startswith("refID: "):
            if transcript:
                tid = transcript["tid"]
                sys.stdout.write("processed transcript {0}.\t\r".format(rid))
                sys.stdout.flush()
                tlist[tid] = transcript.copy()
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
                prof.append([int(pos), int(cnt)])
            transcript["prof"][read_len] = prof
        line = tf.readline()
    if transcript:
        tid = transcript["tid"]
        tlist[tid] = transcript
    tf.close()
    print "\n"
    return tlist

def get_tseq(fn):
    print "getting transcript sequences..."
    ifile = open(fn, 'r')
    tseq = {}
    for rec in SeqIO.parse(ifile, "fasta"):
        tid = rec.id.split()[0]
        seq = str(rec.seq)
        tseq[tid] = seq
    ifile.close()
    return tseq


# read in list of read counts and ref sequence
# output raw read sequences
def get_reads(tlist, tseq, min_len,max_len):
  rseqs = defaultdict(list)
  for tid in tlist:
     for rlen in tlist[tid]["prof"]:
	if rlen < min_len or rlen> max_len:
	  continue
	for pos,count in tlist[tid]["prof"][rlen]:
	  seq = tseq[tid][pos:pos+rlen]
	  rseqs[tid].append((seq,count))
  return rseqs

def get_reads_region(tlist, tseq, min_len,max_len, offset, region_len):
  rseqs = defaultdict(list)
  for tid in tlist:
     for rlen in tlist[tid]["prof"]:
	if rlen < min_len or rlen> max_len:
	  continue
	for pos,count in tlist[tid]["prof"][rlen]:
            if offset < 0:
                seq = tseq[tid][pos+rlen+offset-region_len:pos+rlen+offset]
            else:
                seq = tseq[tid][pos+offset: pos+offset+region_len]
	    rseqs[tid].append((seq,count))
  return rseqs

#=============================
# main
#=============================
def main():
    if len(sys.argv) != 4:
        print "Usage: python extract_read_seq.py input_rlen.hist ref.fasta output.fasta"
        exit(1)
    hist_fn = sys.argv[1]
    rfa = sys.argv[2]
    ofa = sys.argv[3]
    tlist = parse_rlen_hist(hist_fn)
    tseq = get_tseq(rfa)
    reads = get_reads_region(tlist, tseq, 57,61,15,9)
    with open(ofa,'w') as f:
	for tid in reads:
	  for ind,read in enumerate(reads[tid]):
	    for i in range(read[1]):
		f.write('>'+tid+'_'+str(ind)+'_'+str(i)+'\n'+read[0]+'\n')


if __name__ == "__main__": main()
