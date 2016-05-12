#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <seqan/bam_io.h>

using namespace std;

using pos2cnt_t = std::map<size_t, int>;
using rlen2prof_t = std::map<int, pos2cnt_t >;
using rid2rlhist_t = std::map<uint32_t, rlen2prof_t>;
using rid2tlen_t = std::vector<uint32_t>;
using rid2tid_t = std::vector<string>;

void get_ref_len(const seqan::BamHeader& h, rid2tlen_t& rid2tlen, vector<string>& rid2tid)
{
  for (unsigned i = 0; i< length(h.sequenceInfos); ++i) {
    string tid(toCString(h.sequenceInfos[i].i1));
    rid2tlen.emplace_back(h.sequenceInfos[i].i2);
    rid2tid.emplace_back(tid);
  }
}

int get_map_cnt_from_bam_rec(seqan::BamAlignmentRecord& bam_rec)
{
  seqan::BamTagsDict tag(bam_rec.tags);
  unsigned i(0);
  findTagKey(i, tag, seqan::CharString("NH"));
  int val(1);
  extractTagValue(val, tag, i);
  return val;
}

int get_read_len_from_bam_rec(const seqan::BamAlignmentRecord& bam_rec)
{
  int read_len(0);
  bool skip(false);
  for (int i=0; i<length(bam_rec.cigar); ++i) {
    switch (bam_rec.cigar[i].operation) {
    case 'M':
      // add match length to read length
      read_len += bam_rec.cigar[i].count;
      break;
    case 'D':
      // Deletion is a gap to the read
      // read covered one base longer on the transcript
      read_len += 1;
      break;
    case 'N':
      // skipped regions indicate novel splices
      // since reads are mapped to the transcriptome
      // skip these alignments for now
      skip = true;
      break;
    case 'I':
      // Insertion is a gap to the reference
      // no change in read length
    case 'S':
      // softclipping won't change read length
    default:
      break;
    }
    if (skip) break;
  }
  return read_len;
}

int read_len_hist(const char* fn, bool pbegin, rid2tid_t& rid2tid, rid2rlhist_t& rid2hist)
{
  // open bam file
  seqan::BamStream bamIn(fn);
  if (!isGood(bamIn)){
    cerr << "ERROR: Could not open "<<fn<<"!"<<endl;
    exit(1);
  }
  // get reference length info from bam header
  rid2tlen_t rid2tlen;
  get_ref_len(bamIn.header, rid2tlen, rid2tid);
  // get read length and start location info
  seqan::BamAlignmentRecord bam_rec;
  int i=0;
  while(!atEnd(bamIn)){
    if(readRecord(bam_rec,bamIn)!=0){
      cerr << "ERROR: Could not read bam record!\n";
      return true;
    }
    // get mapped reads
    if (!hasFlagUnmapped(bam_rec)){
      // ignore reverse complement alignment
      if (hasFlagRC(bam_rec)) continue;
      // ignore multimappers for now
      if (hasFlagSecondary(bam_rec)) continue;
      // map count in Tag NH, ignore if >1
      if (get_map_cnt_from_bam_rec(bam_rec) > 1) continue;
      // read length is the length on the transcript covered by the read
      int read_len(get_read_len_from_bam_rec(bam_rec));
      // put count in histogram
      if (read_len == 0) continue;
      uint32_t refID(bam_rec.rID);
      int pos_begin(bam_rec.beginPos);
      int pos_end(pos_begin+read_len);
      int pos = pbegin ? pos_begin : pos_end;
      auto it_rid = rid2hist.find(refID);
      // insert refID if not exist
      if (it_rid == rid2hist.end()) {
	rid2hist.emplace(refID, rlen2prof_t{{read_len, pos2cnt_t{{pos, 1}}}});
      }
      else {
	auto it_rlen = it_rid->second.find(read_len);
	// insert read_len if not exist
	if ( it_rlen == it_rid->second.end()) 
	  it_rid->second.emplace(read_len, pos2cnt_t{{pos, 1}});
	// add count if both refID and read_len exist
	else {
	  auto it_pos = it_rlen->second.find(pos);
	  // insert position if not exist
	  if ( it_pos == it_rlen->second.end())
	    it_rlen->second.emplace(pos, 1);
	  else 
	    it_pos->second++;
	}
      }
    }//if(!hasFlagUnmapped(bam_rec))
    if (++i%10000 == 0) {
      cout<<"processed "<<i<<" reads.\r";
      cout.flush();
    }
  }//while(!atEnd(bamIn))
}

void write_rlen_hist(const char* fn, const rid2tid_t& rid2tid, const rid2rlhist_t& rid2hist)
{
  int i = 0;
  ofstream hist_file(fn);
  for (auto it_rid = rid2hist.begin(); it_rid!= rid2hist.end(); ++it_rid) {
    size_t rid(it_rid->first);
    string tid(rid2tid[rid]);
    hist_file<<"refID: "<<rid<<endl;
    hist_file<<"tid: "<<tid<<endl;
    for (auto it_rlen = it_rid->second.begin(); it_rlen != it_rid->second.end(); ++it_rlen)  {
      hist_file<<it_rlen->first<<": ";
      for (auto it_pos = it_rlen->second.begin(); it_pos != it_rlen->second.end(); ++it_pos) {
    	hist_file<<it_pos->first<<" "<<it_pos->second<<", ";
	i += it_pos->second;
      }
      hist_file<<endl;
    }
  }
  hist_file.close();
  cout<<"total transcripts with non-zero coverage: "<<rid2hist.size()<<endl;
  cout<<"total reads: "<<i<<endl;
}

int main(int argc, char ** argv)
{
  if (argc != 3) {
    cout<<"Usage: ./bam_info bam_fname ohist_fname"<<endl;
    exit(1);
  }
  char* bam_in = argv[1];
  char* ofname = argv[2];
  rid2tid_t rid2tid;
  rid2rlhist_t rid2hist;
  cout<<"5'end "<<endl;
  cout<<"parsing bam file..."<<endl;
  read_len_hist(bam_in, true, rid2tid, rid2hist);
  cout<<"writing histogram..."<<endl;
  write_rlen_hist(ofname, rid2tid, rid2hist);
  return 0;
}
