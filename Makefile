
ifeq ($(wildcard /Library),) 

CXX = /opt/rh/devtoolset-2/root/usr/bin/g++
BOOST = /scratch/olson/boost
STDCPP = -std=c++0x 
THREADLIB = -lpthread -lrt

else 

BOOST = /Users/olson/c++/boost
STDCPP = -stdlib=libc++ -std=gnu++11  

endif

default: kser

#OPT = -O2
#OPT = -O2 -pg
OPT = -g
# OPT = -g -DBOOST_ASIO_ENABLE_HANDLER_TRACKING

INC = -I$(BOOST)/include 

CXXFLAGS = $(STDCPP) $(INC) $(OPT)

LDFLAGS  =

LIBS = $(BOOST)/lib/libboost_system.a \
	$(BOOST)/lib/libboost_filesystem.a \
	$(BOOST)/lib/libboost_timer.a \
	$(BOOST)/lib/libboost_chrono.a \
	$(BOOST)/lib/libboost_iostreams.a \
	$(THREADLIB)

x.o: x.cc kguts.h

x: x.o
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

kc.o: kmer.h kserver.h krequest.h klookup.h klookup2.h klookup3.h
klookup2.o: klookup2.h kmer.h global.h
klookup3.o: klookup3.h kmer.h global.h
klookup.o: klookup.h kmer.h global.h
kmer.o: kmer.h
krequest.o: krequest.h kmer.h klookup.h klookup2.h klookup3.h global.h
kser.o: global.h kmer.h kserver.h krequest.h klookup.h klookup2.h klookup3.h
kserver.o: kserver.h kmer.h krequest.h klookup.h klookup2.h klookup3.h
kserver.o: global.h
