CXX = g++
CXXFLAGS = -std=c++11 -O3 -Wno-deprecated-declarations -DSEQAN_HAS_ZLIB=1 -DSEQAN_HAS_BZIP2=1 
INC = -I/home/hw1/.local/include/ -I/opt/local/include
CXXFLAGS += $(INC)
COMPILE.c = $(CXX) $(CXXFLAGS)
LDFLAGS = -lz -lbz2
OUTPUT_OPTION = -o $@

read_len_hist: read_len_hist.o
	$(COMPILE.c) $(OUTPUT_OPTION) $^ $(LDFLAGS)
	chmod u+x $@

collapse_barcode: collapse_barcode.o
	$(COMPILE.c) $(OUTPUT_OPTION) $^ $(LDFLAGS)
	chmod u+x $@

filter_zero_seq: filter_zero_seq.o
	$(COMPILE.c) $(OUTPUT_OPTION) $^ $(LDFLAGS)
	chmod u+x $@
