
ifeq ($(wildcard /Library),) 

#CXX = /opt/rh/devtoolset-2/root/usr/bin/g++
CXX = /disks/olson/gcc-4.9.3/bin/g++
BOOST = /scratch/olson/boost-1.59.0
#BOOST = /scratch/olson/boost
STDCPP = -std=c++14 
THREADLIB = -lpthread -lrt

else 

BOOST = /Users/olson/c++/boost
STDCPP = -stdlib=libc++ -std=gnu++11  

endif

default: kser

# OPT = -O3
OPT = -O2 -pg
OPT = -g -O3
# OPT = -g
# OPT = -g -DBOOST_ASIO_ENABLE_HANDLER_TRACKING

PROFILER_DIR = /scratch/olson/gperftools

#PROFILER_LIB = -L$(PROFILER_DIR)/lib -lprofiler
#PROFILER_INC = -DGPROFILER -I$(PROFILER_DIR)/include

INC = -I$(BOOST)/include 

KMC_DIR = ../KMC/kmc_api
KMC_LIB = $(KMC_DIR)/*.o
KMC_INC = -I$(KMC_DIR)

CXXFLAGS = $(STDCPP) $(INC) $(OPT) $(PROFILER_INC) $(KMC_INC)
CFLAGS = $(INC) $(OPT) $(PROFILER_INC) $(KMC_INC)

LDFLAGS  = -static

LIBS = $(BOOST)/lib/libboost_system.a \
	$(BOOST)/lib/libboost_filesystem.a \
	$(BOOST)/lib/libboost_timer.a \
	$(BOOST)/lib/libboost_chrono.a \
	$(BOOST)/lib/libboost_iostreams.a \
	$(BOOST)/lib/libboost_regex.a \
	$(BOOST)/lib/libboost_thread.a \
	$(BOOST)/lib/libboost_program_options.a \
	$(THREADLIB) \
	$(PROFILER_LIB) \
	$(KMC_LIB)

x.o: x.cc kguts.h

x: x.o
	$(CXX) $(LDFLAGS) $(OPT) -o $@ $^ $(LIBS)

tthr: tthr.o threadpool.o kguts.o
	$(CXX) $(LDFLAGS) $(OPT) -o $@ $^ $(LIBS)

kmerge: kmerge.o
	$(CXX) $(LDFLAGS) $(OPT) -o $@ $^ $(LIBS)

kc: kc.o kmer.o kserver.o krequest.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

kser: kser.o kmer.o kserver.o krequest.o klookup.o klookup2.o klookup3.o kguts.o \
	fasta_parser.o krequest2.o add_request.o threadpool.o matrix_request.o query_request.o lookup_request.o md5.o
	$(CXX) $(LDFLAGS) $(OPT) -o kser $^ $(LIBS)

kfile: kfile.o kguts.o fasta_parser.o
	$(CXX) $(LDFLAGS) $(OPT) -o $@ $^ $(LIBS)

clean:
	rm -f *.o kc kser

depend:
	makedepend -Y *.cc

# DO NOT DELETE

add_request.o: add_request.h compute_request.h krequest2.h kmer.h klookup.h
add_request.o: kguts.h fasta_parser.h klookup2.h klookup3.h threadpool.h
fasta_parser.o: fasta_parser.h
kc.o: kmer.h kserver.h krequest.h klookup.h kguts.h fasta_parser.h klookup2.h
kc.o: klookup3.h krequest2.h compute_request.h threadpool.h
kfile.o: kguts.h fasta_parser.h
kguts.o: kguts.h global.h
klookup2.o: klookup2.h kmer.h kguts.h fasta_parser.h global.h
klookup3.o: klookup3.h kmer.h global.h
klookup.o: klookup.h kmer.h kguts.h fasta_parser.h global.h
kmer.o: popen.h kmer.h
kmerge.o: stringutil.h
krequest2.o: krequest2.h kmer.h klookup.h kguts.h fasta_parser.h klookup2.h
krequest2.o: klookup3.h compute_request.h threadpool.h kserver.h krequest.h
krequest2.o: global.h add_request.h matrix_request.h prot_seq.h
krequest2.o: query_request.h debug.h
krequest.o: krequest.h kmer.h klookup.h kguts.h fasta_parser.h klookup2.h
krequest.o: klookup3.h global.h
kser.o: global.h kmer.h kserver.h krequest.h klookup.h kguts.h fasta_parser.h
kser.o: klookup2.h klookup3.h krequest2.h compute_request.h threadpool.h
kserver.o: kserver.h kmer.h krequest.h klookup.h kguts.h fasta_parser.h
kserver.o: klookup2.h klookup3.h krequest2.h compute_request.h threadpool.h
kserver.o: global.h
lookup_request.o: lookup_request.h compute_request.h krequest2.h kmer.h
lookup_request.o: klookup.h kguts.h fasta_parser.h klookup2.h klookup3.h
lookup_request.o: threadpool.h prot_seq.h global.h
matrix_request.o: matrix_request.h compute_request.h krequest2.h kmer.h
matrix_request.o: klookup.h kguts.h fasta_parser.h klookup2.h klookup3.h
matrix_request.o: threadpool.h prot_seq.h
query_request.o: query_request.h compute_request.h krequest2.h kmer.h
query_request.o: klookup.h kguts.h fasta_parser.h klookup2.h klookup3.h
query_request.o: threadpool.h
threadpool.o: threadpool.h kguts.h
tthr.o: threadpool.h kguts.h
x.o: kguts.h
