#include <unordered_map>
#include <map>
#include <set> 
#include <string> 
#include <cstring>
#include <iostream>
#include <cstdio>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

int main(int argc, char** argv)
{
  if (argc != 4) {
    std::cout<<"Usage: ./collapse_barcode read_in.fq.gz barcode.fq.gz read_out.fq.gz"<<std::endl;
    exit(1);
  }
  char* read_in = argv[1];
  char* barcode = argv[2];
  char* read_out = argv[3];
  seqan::SequenceStream read_stream(read_in), barcode_stream(barcode);
  seqan::SequenceStream out_stream(read_out, seqan::SequenceStream::WRITE);
  // check whether files are opened properly
  std::cout<<read_in<<std::endl;
  if (!isGood(read_stream)) {
    std::cerr << "ERROR: Could not open "<<read_in<<".\n";
    return 1;
  }
  if (!isGood(barcode_stream)) {
    std::cerr << "ERROR: Could not open "<<barcode<<".\n";
    return 1;
  }
  if (!isGood(out_stream)) {
    std::cerr << "ERROR: Could not open "<<read_out<<".\n";
    return 1;
  }

  double tot(0), kept(0), invalid(0), collapse(0);
  std::unordered_map<std::string, std::set<std::string> > barcode2seqs;
  std::map<int, int> barcode_len_hist;
  seqan::CharString read_id, read_seq, qual, barcode_id, barcode_seq;
  // loop through all reads to collapse barcodes
  while (!atEnd(read_stream) and !atEnd(barcode_stream)) {
    // check whether reads are read in properly
    if (readRecord(read_id, read_seq, qual, read_stream) != 0) {
      std::cerr << "ERROR: Could not read from read_in.fq.gz!\n";
      return 1;
    }
    if (readRecord(barcode_id, barcode_seq, barcode_stream) != 0) {
      std::cerr << "ERROR: Could not read from barcode.fq.gz!\n";
      return 1;
    }
    // check whether read id from two files match
    if (strcmp(toCString(read_id), toCString(barcode_id)) != 0) {
      std::cerr<<"ERROR: read id not match barcode id!\n";
      std::cerr<<read_id<<std::endl;
      std::cerr<<barcode_id<<std::endl;
      return 1;
    }
    tot += 1;
    if (int(tot)%10000 == 0) {
      std::cout<<"processed "<<int(tot)<<" reads.\r";
      std::cout.flush();
    }
    std::string barcode(toCString(barcode_seq)), read(toCString(read_seq));
    // check whether read is valid
    if (barcode.length()!=10) {
      auto it = barcode_len_hist.find(barcode.length());
      if (it==barcode_len_hist.end())
	barcode_len_hist.emplace(barcode.length(), 1);
      else
	(it->second)++;
      invalid += 1;
      continue;
    }
    // skip if no read portion before adapter
    if (read.length()==0) {
      invalid += 1;
      continue;
    }
    // check whether barcode and read in record
    // if not, insert record and keep read
    auto it_barcode = barcode2seqs.find(barcode);
    if ( it_barcode != barcode2seqs.end()) {
      auto it_seq = it_barcode->second.find(read);
      // barcode and read in record, skip 
      if ( it_seq != it_barcode->second.end() ) {
	collapse += 1;
	continue;
      }
      else
	it_barcode->second.emplace(read);
    }
    else 
      barcode2seqs.emplace(barcode, std::set<std::string>{read});
    // std::cout << barcode << '\t' << read << '\n';
    if (writeRecord(out_stream, read_id, read_seq, qual) != 0) {
      std::cerr << "ERROR: Could not write to file!\n";
      return 1;
    }
    kept += 1;
  }
  std::printf("total: %.0f kept: %.0f (%.2f %%) collapse: %.0f (%.2f %%) invalid: %.0f (%.2f %%)\n",tot, kept, kept*100/tot, collapse, collapse*100/tot, invalid, invalid*100/tot);
  std::cout<<"invalid barcode histogram: "<<std::endl;
  for (auto i : barcode_len_hist)
    std::cout<<"["<<i.first<<", "<<i.second<<"],";
  return 0;
}
