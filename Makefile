
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
OPT = -g -DBOOST_ASIO_ENABLE_HANDLER_TRACKING

CXXFLAGS = $(STDCPP) -I$(BOOST)/include  $(OPT)

LDFLAGS  =

LIBS = $(BOOST)/lib/libboost_system.a \
	$(BOOST)/lib/libboost_filesystem.a \
	$(BOOST)/lib/libboost_timer.a \
	$(BOOST)/lib/libboost_chrono.a \
	$(THREADLIB)

x: x.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

kc: kc.o kmer.o kserver.o krequest.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

kser: kser.o kmer.o kserver.o krequest.o klookup.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

clean:
	rm *.o kc

# DO NOT DELETE

kc.o: kmer.h kserver.h krequest.h klookup.h
klookup.o: klookup.h kmer.h
kmer.o: /usr/include/stdlib.h /usr/include/Availability.h
kmer.o: /usr/include/AvailabilityInternal.h /usr/include/_types.h
kmer.o: /usr/include/sys/_types.h /usr/include/sys/cdefs.h
kmer.o: /usr/include/sys/_symbol_aliasing.h
kmer.o: /usr/include/sys/_posix_availability.h /usr/include/machine/_types.h
kmer.o: /usr/include/i386/_types.h /usr/include/sys/_pthread/_pthread_types.h
kmer.o: /usr/include/sys/wait.h /usr/include/sys/_types/_pid_t.h
kmer.o: /usr/include/sys/_types/_id_t.h /usr/include/sys/signal.h
kmer.o: /usr/include/sys/appleapiopts.h /usr/include/machine/signal.h
kmer.o: /usr/include/i386/signal.h /usr/include/machine/_mcontext.h
kmer.o: /usr/include/i386/_mcontext.h /usr/include/mach/i386/_structs.h
kmer.o: /usr/include/sys/_pthread/_pthread_attr_t.h
kmer.o: /usr/include/sys/_types/_sigaltstack.h
kmer.o: /usr/include/sys/_types/_ucontext.h
kmer.o: /usr/include/sys/_types/_sigset_t.h /usr/include/sys/_types/_size_t.h
kmer.o: /usr/include/sys/_types/_uid_t.h /usr/include/sys/resource.h
kmer.o: /usr/include/stdint.h /usr/include/sys/_types/_int8_t.h
kmer.o: /usr/include/sys/_types/_int16_t.h /usr/include/sys/_types/_int32_t.h
kmer.o: /usr/include/sys/_types/_int64_t.h /usr/include/_types/_uint8_t.h
kmer.o: /usr/include/_types/_uint16_t.h /usr/include/_types/_uint32_t.h
kmer.o: /usr/include/_types/_uint64_t.h /usr/include/sys/_types/_intptr_t.h
kmer.o: /usr/include/sys/_types/_uintptr_t.h /usr/include/_types/_intmax_t.h
kmer.o: /usr/include/_types/_uintmax_t.h /usr/include/sys/_types/_timeval.h
kmer.o: /usr/include/machine/endian.h /usr/include/i386/endian.h
kmer.o: /usr/include/sys/_endian.h /usr/include/libkern/_OSByteOrder.h
kmer.o: /usr/include/libkern/i386/_OSByteOrder.h /usr/include/alloca.h
kmer.o: /usr/include/sys/_types/_ct_rune_t.h
kmer.o: /usr/include/sys/_types/_rune_t.h /usr/include/sys/_types/_wchar_t.h
kmer.o: /usr/include/sys/_types/_null.h /usr/include/machine/types.h
kmer.o: /usr/include/i386/types.h /usr/include/sys/_types/_dev_t.h
kmer.o: /usr/include/sys/_types/_mode_t.h kmer.h
krequest.o: krequest.h kmer.h klookup.h
kser.o: kmer.h kserver.h krequest.h klookup.h
kserver.o: kserver.h kmer.h krequest.h klookup.h
