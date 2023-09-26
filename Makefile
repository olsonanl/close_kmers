CENTOS=0

KMER_SIZE = 8

ifneq ($(CENTOS),0)
CXX = g++
STDCPP = -std=c++17 -lstdc++fs
THREADLIB = -lpthread -lrt

else
ifeq ($(wildcard /Library),) 

BUILD_TOOLS = /disks/patric-common/runtime/gcc-9.3.0
#BUILD_TOOLS = /disks/patric-common/runtime/gcc-4.9.4
PATH := $(XBUILD_TOOLS)/bin:$(PATH)
export PATH

#CXX = /opt/rh/devtoolset-2/root/usr/bin/g++
#CXX = $(BUILD_TOOLS)/bin/g++
CXX = g++
#BOOST = $(BUILD_TOOLS)
#BOOST = /disks/patric-common/runtime/boost-1.73.0
BOOST = /home/olson/P3/boost/1.83
#STDCPP = -std=c++14 
STDCPP = -std=c++17 -lstdc++fs
THREADLIB = -lpthread -lrt

CXX_LDFLAGS = -Wl,-rpath,$(BUILD_TOOLS)/lib64 -Wl,-rpath,$(BOOST)/lib

else 

BOOST = /Users/olson/c++/boost
STDCPP = -stdlib=libc++ -std=gnu++11  

endif
endif

default: kser

BUILD_DEBUG = 0

ifeq ($(BUILD_DEBUG),1)
OPT = -g -Wall

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

OPT = -O3 -g
#OPT = -O -g -DGLIBCXX_FORCE_NEW
#PROFILE = -pg
#OPT = -O2 $(PROFILE)
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

NUDB = NuDB
NUDB_INCLUDE = -I$(NUDB)/include

#WARNING_FLAGS = -Wconversion -Wall -Werror -Wsign-conversion
CXXFLAGS = -DKMER_SIZE=$(KMER_SIZE) $(STDCPP) $(INC) $(OPT) $(PROFILER_INC) $(KMC_INC) $(BLCR_CFLAGS) $(TBB_CFLAGS) $(NUMA_CFLAGS) $(WARNING_FLAGS) $(LOG4CPP_INC) $(NUDB_INCLUDE)
CFLAGS = $(INC) $(OPT) $(PROFILER_INC) $(KMC_INC) -Wconversion -Wall

# LDFLAGS  = -static
LDFLAGS = $(TBB_LDFLAGS) $(CXX_LDFLAGS) $(PROFILE)

ifeq ($(osrel),CentOS Stream release 8)
BOOST_LIBS = -lboost_system -lboost_log -lboost_filesystem -lboost_timer -lboost_chrono \
	-lboost_iostreams -lboost_regex -lboost_thread -l boost_program_options -lboost_system -llog4cxx
else

#LOG4CPP_DIR = /home/olson/P3/install
LOG4CPP_DIR = /home/olson/P3/install-log4cxx-10
LOG4CPP_LIB = -L$(LOG4CPP_DIR)/lib -llog4cxx -Wl,-rpath=$(LOG4CPP_DIR)/lib
LOG4CPP_INC = -I$(LOG4CPP_DIR)/include

BOOST_LIBS = $(BOOST)/lib/libboost_system.a \
	$(BOOST)/lib/libboost_log.a \
	$(BOOST)/lib/libboost_filesystem.a \
	$(BOOST)/lib/libboost_timer.a \
	$(BOOST)/lib/libboost_chrono.a \
	$(BOOST)/lib/libboost_iostreams.a \
	$(BOOST)/lib/libboost_regex.a \
	$(BOOST)/lib/libboost_thread.a \
	$(BOOST)/lib/libboost_program_options.a \
	$(BOOST)/lib/libboost_system.a 
endif

LIBS = $(BOOST_LIBS) \
	$(THREADLIB) \
	$(PROFILER_LIB) \
	$(BLCR_LIB) \
	$(TBB_LIB) \
	$(NUMA_LIBS) \
	$(LOG4CPP_LIB) \
	-lz

x.o: x.cc kguts.h

x: x.o kmer_encoder.o
	$(CXX) $(LDFLAGS) $(OPT) -o $@ $^ $(LIBS)

x1: x1.o kmer_encoder.o
	$(CXX) $(LDFLAGS) $(OPT) -o $@ $^ $(LIBS)

nudb_test: nudb_test.o nudb_kmer_db.o kmer_nudb.o fasta_parser.o
	$(CXX) $(LDFLAGS) $(OPT) -o $@ $^ $(LIBS)

nudb_call_kmers: nudb_call_kmers.o nudb_kmer_db.o kmer_nudb.o fasta_parser.o
	$(CXX) $(LDFLAGS) $(OPT) -o $@ $^ $(LIBS)

nudb_call_kmers_generic: nudb_call_kmers_generic.o nudb_kmer_db.o fasta_parser.o
	$(CXX) $(LDFLAGS) $(OPT) -o $@ $^ $(LIBS)

nudb_dump: nudb_dump.o
	$(CXX) $(LDFLAGS) $(OPT) -o $@ $^ $(LIBS)

tr: tr.o
	$(CXX) $(LDFLAGS) $(OPT) -o $@ $^ $(LIBS)

load_kmer_nudb: load_kmer_nudb.o kmer_encoder.o
	$(CXX) $(CXX_LDFLAGS) $(LDFLAGS) $(OPT) -o $@ $^ $(LIBS)

xx: xx.o
	$(CXX) $(LDFLAGS) $(OPT) -o $@ $^ $(LIBS)

logtest: logtest.o
	$(CXX) $(LDFLAGS) $(OPT) -o $@ $^ $(LIBS)

tkguts: tkguts.o kguts.o kmer_image.o 
	$(CXX) $(LDFLAGS) $(OPT) -o $@ $^ $(LIBS)

validate_fasta: validate_fasta.o fasta_parser.o
	$(CXX) $(LDFLAGS) $(OPT) -o $@ $^ $(LIBS)

validate_fastq: validate_fastq.o fastq_parser.o trans_table.o dna_seq.o
	$(CXX) $(LDFLAGS) $(OPT) -o $@ $^ $(LIBS)

fastq_to_protein: fastq_to_protein.o fastq_parser.o trans_table.o dna_seq.o
	$(CXX) $(LDFLAGS) $(OPT) -o $@ $^ $(LIBS)

test_family_mapper: test_family_mapper.o fasta_parser.o kguts.o family_mapper.o \
	kmer.o kmer_image.o kmer_inserter.o nr_loader.o threadpool.o kmer_encoder.o
	$(CXX) $(LDFLAGS) $(OPT) -o $@ $^ $(LIBS)

build_signature_kmers: build_signature_kmers.o fasta_parser.o kguts.o kmer_image.o kmer_encoder.o
	$(CXX) $(LDFLAGS) $(OPT) -o $@ $^ $(LIBS)

build_signature_nudb: build_signature_nudb.o fasta_parser.o
	$(CXX) $(LDFLAGS) $(OPT) -o $@ $^ $(LIBS)

recall_proteins: recall_proteins.o fasta_parser.o kmer_nudb.o
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
	fq_process_request.o fastq_parser.o zlib_support.o dna_seq.o trans_table.o family_mapper.o \
	md5.o kmer_image.o numa.o kmer_inserter.o family_reps.o nr_loader.o kmer_encoder.o
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
add_request.o: kguts.h kmer_image.h kmer_types.h kmer_value_types.h
add_request.o: nudb_kmer_db.h kmer_params.h kmer_encoder.h fasta_parser.h
add_request.o: klookup2.h klookup3.h threadpool.h numa.h md5.h
build_signature_kmers.o: nudb_kmer_db.h kmer_value_types.h seed_utils.h
build_signature_kmers.o: operators.h fasta_parser.h kguts.h kmer_image.h
build_signature_kmers.o: kmer_types.h kmer_params.h kmer_encoder.h welford.h
build_signature_kmers.o: global.h function_map.h
build_signature_nudb.o: nudb_kmer_db.h kmer_value_types.h seed_utils.h
build_signature_nudb.o: operators.h fasta_parser.h kmer_types_generic.h
build_signature_nudb.o: welford.h global.h function_map.h
dna_seq.o: dna_seq.h trans_table.h
family_mapper.o: family_mapper.h kguts.h kmer_image.h kmer_types.h
family_mapper.o: kmer_value_types.h nudb_kmer_db.h kmer_params.h
family_mapper.o: kmer_encoder.h kmer.h
family_reps.o: family_reps.h
fasta_parser.o: fasta_parser.h
fastq_parser.o: fastq_parser.h
fastq_to_protein.o: fastq_parser.h trans_table.h dna_seq.h
fq_process_request.o: fq_process_request.h compute_request.h krequest2.h
fq_process_request.o: kmer.h klookup.h kguts.h kmer_image.h kmer_types.h
fq_process_request.o: kmer_value_types.h nudb_kmer_db.h kmer_params.h
fq_process_request.o: kmer_encoder.h fasta_parser.h klookup2.h klookup3.h
fq_process_request.o: threadpool.h numa.h fastq_parser.h prot_seq.h
fq_process_request.o: zlib_support.h dna_seq.h trans_table.h family_mapper.h
fq_process_request.o: kserver.h kmer_inserter.h family_reps.h nr_loader.h
fq_process_request.o: global.h
kc.o: kmer.h kserver.h kmer_inserter.h krequest2.h klookup.h kguts.h
kc.o: kmer_image.h kmer_types.h kmer_value_types.h nudb_kmer_db.h
kc.o: kmer_params.h kmer_encoder.h fasta_parser.h klookup2.h klookup3.h
kc.o: compute_request.h threadpool.h numa.h family_reps.h nr_loader.h
kfile.o: kguts.h kmer_image.h kmer_types.h kmer_value_types.h nudb_kmer_db.h
kfile.o: kmer_params.h kmer_encoder.h fasta_parser.h
kguts.o: kguts.h kmer_image.h kmer_types.h kmer_value_types.h nudb_kmer_db.h
kguts.o: kmer_params.h kmer_encoder.h global.h
klookup.o: klookup.h kmer.h kguts.h kmer_image.h kmer_types.h
klookup.o: kmer_value_types.h nudb_kmer_db.h kmer_params.h kmer_encoder.h
klookup.o: fasta_parser.h global.h
klookup2.o: klookup2.h kmer.h kguts.h kmer_image.h kmer_types.h
klookup2.o: kmer_value_types.h nudb_kmer_db.h kmer_params.h kmer_encoder.h
klookup2.o: fasta_parser.h global.h
klookup3.o: klookup3.h kmer.h global.h
kmer.o: parallel_read.h parallel_read.cc popen.h kmer.h kguts.h kmer_image.h
kmer.o: kmer_types.h kmer_value_types.h nudb_kmer_db.h kmer_params.h
kmer.o: kmer_encoder.h global.h
kmer_encoder.o: kguts.h kmer_image.h kmer_types.h kmer_value_types.h
kmer_encoder.o: nudb_kmer_db.h kmer_params.h kmer_encoder.h
kmer_image.o: kmer_image.h kmer_types.h kmer_value_types.h nudb_kmer_db.h
kmer_image.o: global.h
kmer_inserter.o: kmer_inserter.h kmer.h
kmer_nudb.o: kmer_nudb.h kmer_params.h kmer_types.h kmer_value_types.h
kmer_nudb.o: nudb_kmer_db.h global.h codet.h
kmerge.o: stringutil.h
krequest.o: krequest.h kmer.h klookup.h kguts.h kmer_image.h kmer_types.h
krequest.o: kmer_value_types.h nudb_kmer_db.h kmer_params.h kmer_encoder.h
krequest.o: fasta_parser.h klookup2.h klookup3.h global.h
krequest2.o: krequest2.h kmer.h klookup.h kguts.h kmer_image.h kmer_types.h
krequest2.o: kmer_value_types.h nudb_kmer_db.h kmer_params.h kmer_encoder.h
krequest2.o: fasta_parser.h klookup2.h klookup3.h compute_request.h
krequest2.o: threadpool.h numa.h kserver.h kmer_inserter.h family_reps.h
krequest2.o: nr_loader.h global.h add_request.h matrix_request.h prot_seq.h
krequest2.o: query_request.h lookup_request.h fq_process_request.h
krequest2.o: fastq_parser.h zlib_support.h dna_seq.h trans_table.h
krequest2.o: family_mapper.h debug.h
kser.o: global.h kmer.h kserver.h kmer_inserter.h krequest2.h klookup.h
kser.o: kguts.h kmer_image.h kmer_types.h kmer_value_types.h nudb_kmer_db.h
kser.o: kmer_params.h kmer_encoder.h fasta_parser.h klookup2.h klookup3.h
kser.o: compute_request.h threadpool.h numa.h family_reps.h nr_loader.h
kserver.o: kserver.h kmer.h kmer_inserter.h krequest2.h klookup.h kguts.h
kserver.o: kmer_image.h kmer_types.h kmer_value_types.h nudb_kmer_db.h
kserver.o: kmer_params.h kmer_encoder.h fasta_parser.h klookup2.h klookup3.h
kserver.o: compute_request.h threadpool.h numa.h family_reps.h nr_loader.h
kserver.o: global.h
load_kmer_nudb.o: kmer_encoder.h kmer_params.h kmer_types.h
load_kmer_nudb.o: kmer_value_types.h nudb_kmer_db.h
lookup_request.o: lookup_request.h compute_request.h krequest2.h kmer.h
lookup_request.o: klookup.h kguts.h kmer_image.h kmer_types.h
lookup_request.o: kmer_value_types.h nudb_kmer_db.h kmer_params.h
lookup_request.o: kmer_encoder.h fasta_parser.h klookup2.h klookup3.h
lookup_request.o: threadpool.h numa.h prot_seq.h kserver.h kmer_inserter.h
lookup_request.o: family_reps.h nr_loader.h global.h
matrix_request.o: matrix_request.h compute_request.h krequest2.h kmer.h
matrix_request.o: klookup.h kguts.h kmer_image.h kmer_types.h
matrix_request.o: kmer_value_types.h nudb_kmer_db.h kmer_params.h
matrix_request.o: kmer_encoder.h fasta_parser.h klookup2.h klookup3.h
matrix_request.o: threadpool.h numa.h prot_seq.h
nr_loader.o: nr_loader.h kmer.h threadpool.h kmer_image.h kmer_types.h
nr_loader.o: kmer_value_types.h nudb_kmer_db.h kguts.h kmer_params.h
nr_loader.o: kmer_encoder.h numa.h kmer_inserter.h fasta_parser.h
nudb_call_kmers.o: nudb_kmer_db.h kmer_value_types.h kmer_nudb.h
nudb_call_kmers.o: kmer_params.h kmer_types.h fasta_parser.h
nudb_call_kmers_generic.o: nudb_kmer_db.h kmer_value_types.h kmer_generic.h
nudb_call_kmers_generic.o: kmer_params.h kmer_types_generic.h global.h
nudb_call_kmers_generic.o: codet.h kmer_generic.tcc fasta_parser.h
nudb_dump.o: nudb_kmer_db.h kmer_value_types.h kmer_nudb.h kmer_params.h
nudb_dump.o: kmer_types.h
nudb_kmer_db.o: nudb_kmer_db.h kmer_value_types.h
nudb_test.o: nudb_kmer_db.h kmer_value_types.h kmer_nudb.h kmer_params.h
nudb_test.o: kmer_types.h fasta_parser.h
propagate_names.o: propagate_names.h operators.h
query_request.o: query_request.h compute_request.h krequest2.h kmer.h
query_request.o: klookup.h kguts.h kmer_image.h kmer_types.h
query_request.o: kmer_value_types.h nudb_kmer_db.h kmer_params.h
query_request.o: kmer_encoder.h fasta_parser.h klookup2.h klookup3.h
query_request.o: threadpool.h numa.h
recall_proteins.o: nudb_kmer_db.h kmer_value_types.h kmer_nudb.h
recall_proteins.o: kmer_params.h kmer_types.h seed_utils.h operators.h
recall_proteins.o: fasta_parser.h global.h function_map.h
templ.o: kmer_generic.h kmer_params.h kmer_types_generic.h kmer_value_types.h
templ.o: global.h codet.h kmer_generic.tcc
threadpool.o: threadpool.h kmer_image.h kmer_types.h kmer_value_types.h
threadpool.o: nudb_kmer_db.h kguts.h kmer_params.h kmer_encoder.h numa.h
trans_table.o: trans_table.h
tst_family_reps.o: family_reps.h
unique_prots.o: global.h kguts.h kmer_image.h kmer_types.h kmer_value_types.h
unique_prots.o: nudb_kmer_db.h kmer_params.h kmer_encoder.h fasta_parser.h
validate_fasta.o: fasta_parser.h
validate_fastq.o: fastq_parser.h trans_table.h dna_seq.h
x1.o: kmer_encoder.h kmer_params.h kmer_types.h kmer_value_types.h
x1.o: nudb_kmer_db.h tabsep.h prot_seq.h
zlib_support.o: zlib_support.h
