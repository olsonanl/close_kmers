
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

OPT = -O2
#OPT = -g
# OPT = -g -DBOOST_ASIO_ENABLE_HANDLER_TRACKING

INC = -I$(BOOST)/include 

CXXFLAGS = $(STDCPP) $(INC) $(OPT)

LDFLAGS  =

LIBS = $(BOOST)/lib/libboost_system.a \
	$(BOOST)/lib/libboost_filesystem.a \
	$(BOOST)/lib/libboost_timer.a \
	$(BOOST)/lib/libboost_chrono.a \
	$(THREADLIB)

depend:
	makedepend $(INC) *.cc

x: x.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

kc: kc.o kmer.o kserver.o krequest.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

kser: kser.o kmer.o kserver.o krequest.o klookup.o klookup2.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

clean:
	rm -f *.o kc kser

kc.o: kserver.h krequest.h
klookup.o: kmer.h global.h
krequest.o: kmer.h klookup.h
krequest.o: global.h
kser.o: kserver.h krequest.h klookup.h
kserver.o: global.h
