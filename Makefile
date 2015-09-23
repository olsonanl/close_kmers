
ifeq ($(wildcard /Library),) 

#CXX = /opt/rh/devtoolset-2/root/usr/bin/g++
CXX = /disks/olson/gcc-4.9.3/bin/g++
BOOST = /scratch/olson/boost-1.59.0
#BOOST = /scratch/olson/boost
STDCPP = -std=c++0x 
THREADLIB = -lpthread -lrt

else 

BOOST = /Users/olson/c++/boost
STDCPP = -stdlib=libc++ -std=gnu++11  

endif

default: kser

 OPT = -O3
#OPT = -O2 -pg
# OPT = -g -O2
# OPT = -g
# OPT = -g -DBOOST_ASIO_ENABLE_HANDLER_TRACKING

PROFILER_DIR = /scratch/olson/gperftools

3PROFILER_LIB = -L$(PROFILER_DIR)/lib -lprofiler
#PROFILER_INC = -DGPROFILER -I$(PROFILER_DIR)/include

INC = -I$(BOOST)/include 

KMC_DIR = ../KMC/kmc_api
KMC_LIB = $(KMC_DIR)/*.o
KMC_INC = -I$(KMC_DIR)

CXXFLAGS = $(STDCPP) $(INC) $(OPT) $(PROFILER_INC) $(KMC_INC)

LDFLAGS  = -static

LIBS = $(BOOST)/lib/libboost_system.a \
	$(BOOST)/lib/libboost_filesystem.a \
	$(BOOST)/lib/libboost_timer.a \
	$(BOOST)/lib/libboost_chrono.a \
	$(BOOST)/lib/libboost_iostreams.a \
	$(BOOST)/lib/libboost_regex.a \
	$(BOOST)/lib/libboost_program_options.a \
	$(THREADLIB) \
	$(PROFILER_LIB) \
	$(KMC_LIB)

x.o: x.cc kguts.h

x: x.o
	$(CXX) $(LDFLAGS) $(OPT) -o $@ $^ $(LIBS)

kmerge: kmerge.o
	$(CXX) $(LDFLAGS) $(OPT) -o $@ $^ $(LIBS)

kc: kc.o kmer.o kserver.o krequest.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

kser: kser.o kmer.o kserver.o krequest.o klookup.o klookup2.o klookup3.o kguts.o fasta_parser.o
	$(CXX) $(LDFLAGS) $(OPT) -o kser2 $^ $(LIBS)

kfile: kfile.o kguts.o fasta_parser.o
	$(CXX) $(LDFLAGS) $(OPT) -o $@ $^ $(LIBS)

clean:
	rm -f *.o kc kser

depend:
	makedepend -Y *.cc

# DO NOT DELETE

fasta_parser.o: fasta_parser.h
kc.o: kmer.h kserver.h krequest.h klookup.h kguts.h fasta_parser.h klookup2.h
kc.o: klookup3.h
kfile.o: kguts.h fasta_parser.h
kguts.o: kguts.h
klookup2.o: klookup2.h kmer.h kguts.h fasta_parser.h global.h
klookup3.o: klookup3.h kmer.h global.h
klookup.o: klookup.h kmer.h kguts.h fasta_parser.h global.h
kmer.o: popen.h kmer.h
krequest.o: krequest.h kmer.h klookup.h kguts.h fasta_parser.h klookup2.h
krequest.o: klookup3.h global.h
kser.o: global.h kmer.h kserver.h krequest.h klookup.h kguts.h fasta_parser.h
kser.o: klookup2.h klookup3.h
kserver.o: kserver.h kmer.h krequest.h klookup.h kguts.h fasta_parser.h
kserver.o: klookup2.h klookup3.h global.h
x.o: kguts.h
