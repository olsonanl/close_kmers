
ifeq ($(wildcard /Library),) 

BOOST = /scratch/olson/boost
STDCPP = -std=c++0x 
THREADLIB = -lpthread

else 

BOOST = /Users/olson/c++/boost
STDCPP = -stdlib=libc++ -std=gnu++11  

endif

default: kc

OPT = -O2

CXXFLAGS = $(STDCPP) -I$(BOOST)/include  $(OPT)

LDFLAGS  =


LIBS = $(BOOST)/lib/libboost_system.a $(BOOST)/lib/libboost_filesystem.a $(THREADLIB)

x: x.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

kc: kc.o kmer.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)
