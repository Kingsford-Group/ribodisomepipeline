#include <iostream>
#include <cstdio>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

int main(int argc, char** argv)
{
  if (argc != 3) {
    std::cout<<"Usage: ./filter_zero_seq read_in.fq.gz read_out.fq.gz"<<std::endl;
    exit(1);
  }
  char* read_in = argv[1];
  char* read_out = argv[2];
  seqan::SequenceStream read_stream(read_in);
  seqan::SequenceStream out_stream(read_out, seqan::SequenceStream::WRITE);
  // check whether files are opened properly
  std::cout<<read_in<<std::endl;
  if (!isGood(read_stream)) {
    std::cerr << "ERROR: Could not open "<<read_in<<".\n";
    return 1;
  }
  if (!isGood(out_stream)) {
    std::cerr << "ERROR: Could not open "<<read_out<<".\n";
    return 1;
  }

  int tot(0), kept(0);
  seqan::CharString read_id, read_seq, qual;
  // loop through all reads to collapse barcodes
  while (!atEnd(read_stream)) {
    // check whether reads are read in properly
    if (readRecord(read_id, read_seq, qual, read_stream) != 0) {
      std::cerr << "ERROR: Could not read from read_in.fq.gz!\n";
      return 1;
    }
    if ((++tot)%10000 == 0) {
      std::cout<<"processed "<<tot<<" reads.\r";
      std::cout.flush();
    }
    // skip if no read portion before adapter
    if (length(read_seq)==0)
      continue;
    if (writeRecord(out_stream, read_id, read_seq, qual) != 0) {
      std::cerr << "ERROR: Could not write to file!\n";
      return 1;
    }
    ++kept;
  }
  std::printf("total: %d kept: %d (%.2f %%).\n", tot, kept, kept*100.0/tot);
  return 0;
}
