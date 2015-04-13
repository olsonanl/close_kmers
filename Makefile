BOOST = /scratch/olson/boost
#BOOST = /Users/olson/c++/boost

default: kc
OPT = -O2

STDCPP = -std=c++0x 
#STDCPP = -stdlib=libc++ -std=gnu++11  

CXXFLAGS = $(STDCPP) -I$(BOOST)/include  $(OPT)

LDFLAGS  = -L/Users/olson/c++/boost/lib 

THREADLIB = -lpthread

LIBS = $(BOOST)/lib/libboost_system.a $(BOOST)/lib/libboost_filesystem.a $(THREADLIB)

x: x.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

kc: kc.o kmer.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)
