
ifeq ($(wildcard /Library),) 

BUILD_TOOLS = /disks/patric-common/runtime/gcc-4.9.4
PATH := $(BUILD_TOOLS)/bin:$(PATH)
export PATH

#CXX = /opt/rh/devtoolset-2/root/usr/bin/g++
CXX = $(BUILD_TOOLS)/bin/g++
BOOST = $(BUILD_TOOLS)
#BOOST = /scratch/olson/boost
STDCPP = -std=c++14 
THREADLIB = -lpthread -lrt

CXX_LDFLAGS = -Wl,-rpath,$(BUILD_TOOLS)/lib64

else 

BOOST = /Users/olson/c++/boost
STDCPP = -stdlib=libc++ -std=gnu++11  

endif

default: kser

BUILD_DEBUG = 1

ifeq ($(BUILD_DEBUG),1)
OPT = -g

# OPT = -g -DBOOST_ASIO_ENABLE_HANDLER_TRACKING

else
ifeq ($(BUILD_DEBUG),2)
UFLAGS = \
-u malloc \
-u free \
-u realloc \
-u getenv \
-u setenv \
-u __errno_location \
-u pthread_key_create \
-u pthread_key_delete \
-u pthread_setspecific \
-u pthread_getspecific \
-u pthread_spin_init \
-u pthread_spin_destroy \
-u pthread_spin_lock \
-u pthread_spin_trylock \
-u pthread_spin_unlock \
-u pthread_mutex_init \
-u pthread_mutex_destroy \
-u pthread_mutex_trylock \
-u pthread_mutex_lock \
-u pthread_mutex_unlock \
-u pthread_cond_init \
-u pthread_cond_destroy \
-u pthread_cond_signal \
-u pthread_cond_wait \
-u _pthread_cleanup_push \
-u _pthread_cleanup_pop \
-u pthread_setcancelstate \
-u pthread_self \
-u pthread_yield


OPT = -g -O3 -DTBB_USE_THREADING_TOOLS -fno-omit-frame-pointer $(UFLAGS)


# OPT = -g -DBOOST_ASIO_ENABLE_HANDLER_TRACKING

else

OPT = -O3
#OPT = -O2 -pg
#OPT = -g -O3

endif
endif

PROFILER_DIR = /scratch/olson/gperftools

#PROFILER_LIB = -L$(PROFILER_DIR)/lib -lprofiler
#PROFILER_INC = -DGPROFILER -I$(PROFILER_DIR)/include

INC = -I$(BOOST)/include 

KMC_DIR = ../KMC/kmc_api
KMC_LIB = $(KMC_DIR)/*.o
KMC_INC = -I$(KMC_DIR)

LOG4CPP_DIR = /home/olson/P3/install
LOG4CPP_LIB = -L$(LOG4CPP_DIR)/lib -llog4cxx -Wl,-rpath=$(LOG4CPP_DIR)/lib
LOG4CPP_INC = -I$(LOG4CPP_DIR)/include

USE_TBB = 1
USE_NUMA = 0

ifneq ($(USE_TBB),0)

TBB_DEBUG = 0
TBB_INC_DIR = $(BUILD_TOOLS)/include

ifeq ($(TBB_DEBUG),1)
TBB_LIB_DIR = $(BUILD_TOOLS)/lib
#TBB_LIB_DIR = /disks/olson/checkpoint/tbb44_20150728oss/build/linux_intel64_gcc_cc4.9.3_libc2.12_kernel2.6.32_debug
TBB_LIBS = -ltbbmalloc_debug -ltbb_debug
else
TBB_LIB_DIR = $(BUILD_TOOLS)/lib
#TBB_LIB_DIR = /disks/olson/checkpoint/tbb44_20150728oss/build/linux_intel64_gcc_cc4.9.3_libc2.12_kernel2.6.32_release
TBB_LIBS = -ltbbmalloc -ltbb
endif

ifneq ($(TBB_INC_DIR),"")
TBB_CFLAGS = -DUSE_TBB -I$(TBB_INC_DIR)
TBB_LIB = -L$(TBB_LIB_DIR) $(TBB_LIBS)
TBB_LDFLAGS = -Wl,-rpath,$(TBB_LIB_DIR) 
endif

endif

ifneq ($(USE_NUMA),0)

NUMA_LIBS = -lhwloc
NUMA_CFLAGS = -DUSE_NUMA

endif

#BLCR_DIR = /disks/patric-common/blcr

ifneq ($(BLCR_DIR),)
BLCR_CFLAGS = -DBLCR_SUPPORT -I$(BLCR_DIR)/include
BLCR_LIB = -L$(BLCR_DIR)/lib -lcr
endif

#WARNING_FLAGS = -Wconversion -Wall -Werror
CXXFLAGS = $(STDCPP) $(INC) $(OPT) $(PROFILER_INC) $(KMC_INC) $(BLCR_CFLAGS) $(TBB_CFLAGS) $(NUMA_CFLAGS) $(WARNING_FLAGS) $(LOG4CPP_INC)
CFLAGS = $(INC) $(OPT) $(PROFILER_INC) $(KMC_INC) -Wconversion -Wall

# LDFLAGS  = -static
LDFLAGS = $(TBB_LDFLAGS) $(CXX_LDFLAGS)

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
	$(BLCR_LIB) \
	$(TBB_LIB) \
	$(NUMA_LIBS) \
	$(LOG4CPP_LIB)

x.o: x.cc kguts.h

x: x.o
	$(CXX) $(LDFLAGS) $(OPT) -o $@ $^ $(LIBS)

xx: xx.o
	$(CXX) $(LDFLAGS) $(OPT) -o $@ $^ $(LIBS)

tkguts: tkguts.o kguts.o kmer_image.o 
	$(CXX) $(LDFLAGS) $(OPT) -o $@ $^ $(LIBS)

validate_fasta: validate_fasta.o fasta_parser.o
	$(CXX) $(LDFLAGS) $(OPT) -o $@ $^ $(LIBS)

validate_fastq: validate_fastq.o fastq_parser.o trans_table.o dna_seq.o
	$(CXX) $(LDFLAGS) $(OPT) -o $@ $^ $(LIBS)

tt: tt.o
	$(CXX) $(LDFLAGS) $(OPT) -o $@ $^ $(LIBS)

mm: mm.o
	$(CXX) $(LDFLAGS) $(OPT) -o $@ $^ $(LIBS)

tthr: tthr.o threadpool.o kguts.o
	$(CXX) $(LDFLAGS) $(OPT) -o $@ $^ $(LIBS)

kmerge: kmerge.o
	$(CXX) $(LDFLAGS) $(OPT) -o $@ $^ $(LIBS) $(KMC_LIB)

kc: kc.o kmer.o kserver.o krequest.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

kser: kser.o kmer.o kserver.o kguts.o \
	fasta_parser.o krequest2.o add_request.o threadpool.o matrix_request.o query_request.o lookup_request.o \
	md5.o kmer_image.o numa.o kmer_inserter.o family_reps.o
	$(CXX) $(LDFLAGS) $(OPT) -o kser $^ $(LIBS)

tst_family_reps: family_reps.o tst_family_reps.o
	$(CXX) $(LDFLAGS) $(OPT) -o $@ $^ $(LIBS)

unique_prots: unique_prots.o kguts.o kmer_image.o fasta_parser.o
	$(CXX) $(LDFLAGS) $(OPT) -o $@ $^ $(LIBS)

kfile: kfile.o kguts.o fasta_parser.o
	$(CXX) $(LDFLAGS) $(OPT) -o $@ $^ $(LIBS)

propagate_names: propagate_names.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)
clean:
	rm -f *.o kc kser

depend:
	makedepend -Y *.cc

# DO NOT DELETE

add_request.o: add_request.h compute_request.h krequest2.h kmer.h klookup.h
add_request.o: kguts.h kmer_image.h fasta_parser.h klookup2.h klookup3.h
add_request.o: threadpool.h numa.h md5.h
family_reps.o: family_reps.h
fasta_parser.o: fasta_parser.h
fastq_parser.o: fastq_parser.h
kc.o: kmer.h kserver.h kmer_inserter.h krequest2.h klookup.h kguts.h
kc.o: kmer_image.h fasta_parser.h klookup2.h klookup3.h compute_request.h
kc.o: threadpool.h numa.h family_reps.h
kfile.o: kguts.h kmer_image.h fasta_parser.h
kguts.o: kguts.h kmer_image.h global.h
klookup2.o: klookup2.h kmer.h kguts.h kmer_image.h fasta_parser.h global.h
klookup3.o: klookup3.h kmer.h global.h
klookup.o: klookup.h kmer.h kguts.h kmer_image.h fasta_parser.h global.h
kmer.o: parallel_read.h parallel_read.cc popen.h kmer.h kguts.h kmer_image.h
kmer.o: global.h
kmerge.o: stringutil.h
kmer_image.o: kmer_image.h global.h
kmer_inserter.o: kmer_inserter.h kmer.h
krequest2.o: krequest2.h kmer.h klookup.h kguts.h kmer_image.h fasta_parser.h
krequest2.o: klookup2.h klookup3.h compute_request.h threadpool.h numa.h
krequest2.o: kserver.h kmer_inserter.h family_reps.h global.h add_request.h
krequest2.o: matrix_request.h prot_seq.h query_request.h lookup_request.h
krequest2.o: debug.h
krequest.o: krequest.h kmer.h klookup.h kguts.h kmer_image.h fasta_parser.h
krequest.o: klookup2.h klookup3.h global.h
kser.o: global.h kmer.h kserver.h kmer_inserter.h krequest2.h klookup.h
kser.o: kguts.h kmer_image.h fasta_parser.h klookup2.h klookup3.h
kser.o: compute_request.h threadpool.h numa.h family_reps.h
kserver.o: kserver.h kmer.h kmer_inserter.h krequest2.h klookup.h kguts.h
kserver.o: kmer_image.h fasta_parser.h klookup2.h klookup3.h
kserver.o: compute_request.h threadpool.h numa.h family_reps.h global.h
lookup_request.o: lookup_request.h compute_request.h krequest2.h kmer.h
lookup_request.o: klookup.h kguts.h kmer_image.h fasta_parser.h klookup2.h
lookup_request.o: klookup3.h threadpool.h numa.h prot_seq.h kserver.h
lookup_request.o: kmer_inserter.h family_reps.h global.h
matrix_request.o: matrix_request.h compute_request.h krequest2.h kmer.h
matrix_request.o: klookup.h kguts.h kmer_image.h fasta_parser.h klookup2.h
matrix_request.o: klookup3.h threadpool.h numa.h prot_seq.h
propagate_names.o: propagate_names.h operators.h
query_request.o: query_request.h compute_request.h krequest2.h kmer.h
query_request.o: klookup.h kguts.h kmer_image.h fasta_parser.h klookup2.h
query_request.o: klookup3.h threadpool.h numa.h
threadpool.o: threadpool.h kmer_image.h kguts.h numa.h
trans_table.o: trans_table.h
tst_family_reps.o: family_reps.h
unique_prots.o: global.h kguts.h kmer_image.h fasta_parser.h
validate_fasta.o: fasta_parser.h
validate_fastq.o: fastq_parser.h trans_table.h dna_seq.h
x.o: kguts.h kmer_image.h
