#include "kguts.h"

#include <stdio.h>
#include <stdlib.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <ctype.h>
#include <fcntl.h>
#include <unistd.h>
#include <signal.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <errno.h>
#include <string.h>

static const  char genetic_code[64] = {
    'K','N','K','N','T','T','T','T','R','S','R','S','I','I','M','I',
    'Q','H','Q','H','P','P','P','P','R','R','R','R','L','L','L','L',
    'E','D','E','D','A','A','A','A','G','G','G','G','V','V','V','V',
    '*','Y','*','Y','S','S','S','S','*','C','W','C','L','F','L','F'
};

static const  char prot_alpha[20] = {
    'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y' 
};

KmerGuts::KmerGuts(const std::string &data_dir)
{
    tot_lookups = 0;
    retry  = 0;

    debug = 0;
    aa = 0;
    hits_only = 0;
    size_hash = 1400303159;
    write_mem_map = 0;

    num_hits = 0;
    num_oI = 0;

    order_constraint = 0;
    min_hits = 5;
    min_weighted_hits = 0;
    max_gap  = 200;

    pIseq = (unsigned char *) malloc(MAX_SEQ_LEN / 3);
    cdata = (char *) malloc(MAX_SEQ_LEN);
    data = (char *) malloc(MAX_SEQ_LEN);
    pseq = (char *) malloc(MAX_SEQ_LEN / 3);

    kmersH = init_kmers(data_dir.c_str());
}


/* =========================== end of reduction global variables ================= */

unsigned char KmerGuts::to_amino_acid_off(char c) {
  switch (c)
    {
      case 'A':
	return 0;

      case 'C':
	return 1;

      case 'D':
	return 2;

      case 'E':
	return 3;

      case 'F':
	return 4;

      case 'G':
	return 5;

      case 'H':
	return 6;

      case 'I':
	return 7;

      case 'K':
	return 8;

      case 'L':
	return 9;

      case 'M':
	return 10;

      case 'N':
	return 11;

      case 'P':
	return 12;

      case 'Q':
	return 13;

      case 'R':
	return 14;

      case 'S':
	return 15;

      case 'T':
	return 16;

      case 'V':
	return 17;

      case 'W':
	return 18;

      case 'Y':
	return 19;

      default:
	return 20;
    }
}

char KmerGuts::comp(char c)
{
    switch (c)
    {
      case 'a':
	return 't';
      case 'A':
	return 'T';

      case 'c':
	return 'g';
      case 'C':
	return 'G';

      case 'g':
	return 'c';
      case 'G':
	return 'C';

      case 't':
      case 'u':
	return 'a';
      case 'T':
      case 'U':
	return 'A';

      case 'm':
	return 'k';
      case 'M':
	return 'K';

      case 'r':
	return 'y';
      case 'R':
	return 'Y';

      case 'w':
	return 'w';
      case 'W':
	return 'W';

      case 's':
	return 'S';
      case 'S':
	return 'S';

      case 'y':
	return 'r';
      case 'Y':
	return 'R';

      case 'k':
	return 'm';
      case 'K':
	return 'M';

      case 'b':
	return 'v';
      case 'B':
	return 'V';

      case 'd':
	return 'h';
      case 'D':
	return 'H';

      case 'h':
	return 'd';
      case 'H':
	return 'D';

      case 'v':
	return 'b';
      case 'V':
	return 'B';

      case 'n':
	return 'n';
      case 'N':
	return 'N';

      default:
	return c;
    }
}

void KmerGuts::rev_comp(const char *data,char *cdata) {

  int n = strlen(data);
  const char *p  = data + (n-1);
  char *pc = cdata;
  while (n--) {
    *(pc++) = compl(*(p--));
  }
  *pc = 0;
}

unsigned long long KmerGuts::encoded_kmer(unsigned char *p) {
  unsigned long long encodedK = *p;
  int i;
  for (i=1; (i <= K-1); i++) {
    encodedK = (encodedK * 20) + *(p+i);
  }

  if (encodedK > MAX_ENCODED) {
    fprintf(stderr,"bad encoding - input must have included invalid characters\n");
    for (i=0; (i < K); i++) {
      fprintf(stderr,"%d ",*(p+i));
    }
    fprintf(stderr,"\n");
    exit(2);
  }
  return encodedK;
}

unsigned long long KmerGuts::encoded_aa_kmer(char *p)
{
  unsigned char aa_off[K];
  int j;
  for (j=0; (j < K); j++) {
    int prot_c = *(p+j);
    aa_off[j] = to_amino_acid_off(prot_c);
  }
  return encoded_kmer(aa_off);
}

void KmerGuts::decoded_kmer(unsigned long long encodedK,char *decoded) {
  
  int i;
  *(decoded+K) = '\0';
  unsigned long long x = encodedK;

  for (i=K-1; (i >= 0); i--) {
    *(decoded+i) = prot_alpha[x % 20];
    x = x / 20;
  }
}


int KmerGuts::dna_char(char c)
{
    switch (c)
    {
      case 'a':
      case 'A':
	return 0;

      case 'c':
      case 'C':
	return 1;

      case 'g':
      case 'G':
	return 2;

      case 't':
      case 'u':
      case 'T':
      case 'U':
	return 3;

      default:
	return 4;;
    }
}

void KmerGuts::translate(const char *seq,int off,char *pseq, unsigned char *pIseq) {

  int i;
  int max = strlen(seq) - 3;
  char *p = pseq;
  unsigned char *pI = pIseq;
  for (i=off; (i <= max); ) {
    int c1 = dna_char(seq[i++]);
    int c2 = dna_char(seq[i++]);
    int c3 = dna_char(seq[i++]);
    if ((c1 < 4) && (c2 < 4) && (c3 < 4)) {
      int I = (c1 * 16) + (c2 * 4) + c3;
      char prot_c = genetic_code[I];
      *(p++) = prot_c;
      *(pI++) = to_amino_acid_off(prot_c);
    }
    else {
      *(p++)  = 'x';
      *(pI++) = 20;
    }
  }
  *p = 0;
  *pI = 21;
  if (debug >= 3) {
    fprintf(stderr,"len-seq=%d max=%d p=%ld\n",(int) strlen(seq),max,p-pseq);
  }
}

#define MAX_FUNC_OI_INDEX 1000000
#define MAX_FUNC_OI_VALS  100000000

char **KmerGuts::load_indexed_ar(char *filename,int *sz) {
    char **index_ar = (char **) malloc(MAX_FUNC_OI_INDEX * sizeof(char *));
    char *vals      = (char *)malloc(MAX_FUNC_OI_VALS);
  char *p         = vals;
  FILE *ifp      = fopen(filename,"r");
  if (ifp == NULL) { 
    fprintf(stderr,"could not open %s\n",filename);
    exit(1);
  }

  *sz = 0;
  int j;
  while ((fscanf(ifp,"%d\t",&j) == 1) && fgets(p,1000,ifp)) {
    if (*sz != j) {
      fprintf(stderr,"Your index must be dense and in order (see line %ld, should be %d)\n",p-vals,*sz);
      exit(1);
    }
    /* fprintf(stderr,"%d is %s\n",*sz,p); */
    index_ar[*sz] = p;               /* the fgets leaves the \n at the end of each line */
    p += strlen(index_ar[*sz]) -1;
    *(p++) = '\0';
    if ((*sz >= MAX_FUNC_OI_INDEX) || ((p-vals) > (MAX_FUNC_OI_VALS - 1000))) {
      fprintf(stderr,"Your function or oI index arrays are too small; bump MAX_FUNC_OI_INDEX and MAX_FUNC_OI_VALS\n");
      exit(1);
    }

    *sz += 1;
  }
  return index_ar;
}

char **KmerGuts::load_functions(char *file) {
  int sz;
  return load_indexed_ar(file,&sz);
}

char **KmerGuts::load_otus(char *file) {
  int sz;
  return load_indexed_ar(file,&sz);
}

long long KmerGuts::find_empty_hash_entry(sig_kmer_t sig_kmers[],unsigned long long encodedK) {
    long long hash_entry = encodedK % size_hash;
    while (sig_kmers[hash_entry].which_kmer <= MAX_ENCODED)
      hash_entry = (hash_entry+1)%size_hash;
    return hash_entry;
}

long long KmerGuts::lookup_hash_entry(sig_kmer_t sig_kmers[],unsigned long long encodedK) {
    long long  hash_entry = encodedK % size_hash;
    if (debug >= 2)
      tot_lookups++;
    while ((sig_kmers[hash_entry].which_kmer <= MAX_ENCODED) && (sig_kmers[hash_entry].which_kmer != encodedK)) {
      if (debug >= 2)
	retry++;
      hash_entry++;
      if (hash_entry == size_hash)
	hash_entry = 0;
    }
    if (sig_kmers[hash_entry].which_kmer > MAX_ENCODED) {
      return -1;
    }
    else {
      return hash_entry;
    }
}

KmerGuts::kmer_memory_image_t *KmerGuts::load_raw_kmers(char *file,unsigned long long num_entries, unsigned long long *alloc_sz) {
  /*
   * Allocate enough memory to hold the kmer_memory_image_t header plus the hash table itself.
   */
  *alloc_sz = sizeof(kmer_memory_image_t) + (sizeof(sig_kmer_t) * num_entries);

  kmer_memory_image_t *image = (kmer_memory_image_t *) malloc(*alloc_sz);

  /*
   * Initialize our table pointer to the first byte after the header.
   */
  image->num_sigs = num_entries;
  image->entry_size = sizeof(sig_kmer_t);
  image->version = (long long) VERSION;

  sig_kmer_t *sig_kmers = (sig_kmer_t *) (image + 1);

  FILE *ifp      = fopen(file,"r");
  if (ifp == NULL) { 
    fprintf(stderr,"could not open %s",file);
    exit(1);
  }

  long long i;
  for (i=0; (i < size_hash); i++)
    sig_kmers[i].which_kmer = MAX_ENCODED + 1;

  char kmer_string[K+1];
  int end_off;
  int fI;
  float f_wt;
  int oI;
  long long loaded = 0;
  while (fscanf(ifp,"%s\t%d\t%d\t%f\t%d",
		kmer_string,&end_off,&fI,&f_wt,&oI) >= 4) {
    unsigned long long encodedK = encoded_aa_kmer(kmer_string);
    long long hash_entry = find_empty_hash_entry(sig_kmers,encodedK);
    loaded++;
    if (loaded >= (size_hash / 2)) {
      fprintf(stderr,"Your Kmer hash is half-full; use -s (and -w) to bump it\n");
      exit(1);
    }
    sig_kmers[hash_entry].which_kmer     = encodedK;
    sig_kmers[hash_entry].avg_from_end   = end_off;
    sig_kmers[hash_entry].function_index = fI;
    sig_kmers[hash_entry].otu_index      = oI;
    sig_kmers[hash_entry].function_wt    = f_wt;
  }
  if (debug >= 2)
    fprintf(stderr,"loaded %lld kmers\n",loaded);

  return image;
}

KmerGuts::kmer_handle_t *KmerGuts::init_kmers(const char *dataD) {
    kmer_handle_t *handle = (kmer_handle_t *) malloc(sizeof(kmer_handle_t));

  kmer_memory_image_t *image;

  char file[300];
  strcpy(file,dataD);
  strcat(file,"/function.index");
  handle->function_array = load_functions(file);

  strcpy(file,dataD);
  strcat(file,"/otu.index");
  handle->otu_array      = load_otus(file);

  char fileM[300];
  strcpy(fileM,dataD);
  strcat(fileM,"/kmer.table.mem_map");

  if (write_mem_map) {
    unsigned long long sz, table_size;
    strcpy(file,dataD);
    strcat(file,"/final.kmers");
    
    unsigned long long image_size;

    image = load_raw_kmers(file, size_hash, &image_size);

    handle->kmer_table = (sig_kmer_t *) (image + 1);
    handle->num_sigs   = image->num_sigs;

    FILE *fp = fopen(fileM,"w");
    if (fp == NULL) { 
      fprintf(stderr,"could not open %s for writing: %s ",fileM, strerror(errno));
      exit(1);
    }
    fwrite(image, image_size, 1, fp);
    fclose(fp);

    strcpy(fileM,dataD);
    strcat(fileM,"/size_hash.and.table_size");
    fp = fopen(fileM,"w");
    fprintf(fp,"%lld\t%lld\n",sz,table_size);
    fclose(fp);
  }
  else {
    int fd;
    if ((fd = open(fileM, O_RDONLY)) == -1) {
      perror("open");
      exit(1);
    }

    /*
     * Set up for creating memory image from file. Start by determining file size
     * on disk with a stat() call.
     */
    struct stat sbuf;
    if (stat(fileM, &sbuf) == -1) {
      fprintf(stderr, "stat %s failed: %s\n", fileM, strerror(errno));
      exit(1);
    }
    unsigned long long file_size = sbuf.st_size;

    /* 
     * Memory map.
     */
    int flags = MAP_SHARED;
    #ifdef MAP_POPULATE
    flags |= MAP_POPULATE;
    #endif
    
    image = (kmer_memory_image_t *) mmap((caddr_t)0, file_size, PROT_READ, flags, fd, 0);

    if (image == (kmer_memory_image_t *)(-1)) {
      fprintf(stderr, "mmap of kmer_table %s failed: %s\n", fileM, strerror(errno));
      exit(1);
    }

    /* 
     * Our image is mapped. Validate against the current version of this code.
     */
    if (image->version != (long long) VERSION) {
      fprintf(stderr, "Version mismatch for file %s: file has %lld code has %lld\n", 
	      fileM, image->version, (long long) VERSION);
      exit(1);
    }

    if (image->entry_size != (unsigned long long) sizeof(sig_kmer_t)) {
      fprintf(stderr, "Version mismatch for file %s: file has entry size %lld code has %lld\n",
	      fileM, image->entry_size, (unsigned long long) sizeof(sig_kmer_t));
      exit(1);
    }

    size_hash = image->num_sigs;
    handle->num_sigs = size_hash;
    handle->kmer_table = (sig_kmer_t *) (image + 1);

    /* Validate overall file size vs the entry size and number of entries */
    if (file_size != ((sizeof(sig_kmer_t) * image->num_sigs) + sizeof(kmer_memory_image_t))) {
      fprintf(stderr, "Version mismatch for file %s: file size does not match\n", fileM);
      exit(1);
    }

    fprintf(stderr, "Set size_hash=%lld from file size %lld\n", size_hash, file_size);

  }
  return handle;
}

void KmerGuts::advance_past_ambig(unsigned char **p,unsigned char *bound) {

  if (K == 5) {
    while (((*p) < bound) &&
	   ((*(*p) == 20)     || 
            (*((*p)+1) == 20) || 
            (*((*p)+2) == 20) || 
            (*((*p)+3) == 20) || 
	    (*((*p)+4) == 20) )) {
      (*p)++;
    }
  }
  else {   /*  ##### ASSUMING K == 8 #### */
    int bad = 1;
    while ((*p < bound) && (bad == 1)) {
      bad = 0;
      if      (*((*p)+7) == 20) {
	bad = 1;
	(*p) += 8;
      }
      else if (*((*p)+6) == 20) {
	bad = 1;
	(*p) += 7;
      }
      else if (*((*p)+5) == 20) {
	bad = 1;
	(*p) += 6;
      }
      else if (*((*p)+4) == 20) {
	bad = 1;
	(*p) += 5;
      }
      else if (*((*p)+3) == 20) {
	bad = 1;
	(*p) += 4;
      }
      else if (*((*p)+2) == 20) {
	bad = 1;
	(*p) += 3;
      }
      else if (*((*p)+1) == 20) {
	bad = 1;
	(*p) += 2;
      }
      else if (*((*p)+0) == 20) {
	bad = 1;
	(*p) += 1;
      }
    } 
  }
}

void KmerGuts::process_set_of_hits(std::shared_ptr<std::vector<KmerCall>> calls,
				   std::shared_ptr<KmerOtuStats> otu_stats)
{
    if (!otu_stats && !calls)
	return;
    
    int fI_count = 0;
    float weighted_hits = 0;
    int last_hit=0;
    int i=0;
    while (i < num_hits) {
	if (hits[i].fI == current_fI) {
	    last_hit = i;
	    fI_count++;
	    weighted_hits += hits[i].function_wt;
	}
	i++;
    }
    if ((fI_count >= min_hits) && (weighted_hits >= min_weighted_hits))
    {
	if (calls)
	    calls->push_back({ hits[0].from0_in_prot, hits[last_hit].from0_in_prot+(K-1), fI_count, current_fI, weighted_hits });

	/* once we have decided to call a region, we take the kmers for fI and
	   add them to the counts maintained to assign an OTU to the sequence */

	if (otu_stats)
	{
	    for (i=0; (i <= last_hit); i++)
	    {
		if (hits[i].fI == current_fI)
		{
		    otu_stats->otu_map[hits[i].oI]++;
		}
	    }
	}
    }

    if ((hits[num_hits-2].fI != current_fI) && (hits[num_hits-2].fI == hits[num_hits-1].fI)) {
	current_fI = hits[num_hits-1].fI;
	hits[0] = hits[num_hits - 2];
	hits[1] = hits[num_hits - 1];
	num_hits                 = 2;
    }
    else {
	num_hits = 0;
    }
}

void KmerGuts::gather_hits(int ln_DNA, char strand,int prot_off,const char *pseq,
			   unsigned char *pIseq,
			   std::shared_ptr<std::vector<KmerCall>> calls,
			   std::function<void(sig_kmer_t &)> hit_cb,
			   std::shared_ptr<KmerOtuStats> otu_stats)
{
    unsigned char *p = pIseq;
    /* pseq and pIseq are the same length */

    unsigned char *bound = pIseq + strlen(pseq) - K;
    advance_past_ambig(&p,bound);
    unsigned long long encodedK=0;
    if (p < bound) {
	encodedK = encoded_kmer(p);
    }
    while (p < bound) {
	long long  where = lookup_hash_entry(kmersH->kmer_table,encodedK);
	if (where >= 0) {
	    sig_kmer_t *kmers_hash_entry = &(kmersH->kmer_table[where]);
	    int avg_off_end = kmers_hash_entry->avg_from_end;
	    int fI        = kmers_hash_entry->function_index;
	    int oI          = kmers_hash_entry->otu_index;
	    float f_wt      = kmers_hash_entry->function_wt;
	    if (hit_cb)
		hit_cb(*kmers_hash_entry);

	    if ((num_hits > 0) && (hits[num_hits-1].from0_in_prot + max_gap) < (p-pIseq))
	    {
		if (num_hits >= min_hits)
		{
		    process_set_of_hits(calls, otu_stats);
		}
		else
		{
		    num_hits = 0;
		}
	    }

	    if (num_hits == 0)
	    {
		current_fI = fI;   /* if this is the first, set the current_fI */
	    }

	    if ((! order_constraint) || (num_hits == 0) ||
		((fI == hits[num_hits-1].fI) &&
		 (abs(((p-pIseq) - hits[num_hits-1].from0_in_prot) -
		      (hits[num_hits-1].avg_off_from_end - avg_off_end)
		     ) <= 20)))
	    {
		/* we have a new hit, so we add it to the global set of hits */
		hits[num_hits].oI = oI;
		hits[num_hits].fI = fI;
		hits[num_hits].from0_in_prot = p-pIseq;
		hits[num_hits].avg_off_from_end = avg_off_end;
		hits[num_hits].function_wt = f_wt;
		if (num_hits < MAX_HITS_PER_SEQ - 2)
		    num_hits++;
		if ((num_hits > 1) && (current_fI != fI) &&           /* if we have a pair of new fIs, it is time to */
		    (hits[num_hits-2].fI == hits[num_hits-1].fI))    /* process one set and initialize the next */
		{
		    process_set_of_hits(calls, otu_stats);
		}
	    }
	}
	p++;
	if (p < bound) {
	    if (*(p+K-1) < 20) {
		encodedK = ((encodedK % CORE) * 20L) + *(p+K-1);
	    }
	    else {
		p += K;
		advance_past_ambig(&p,bound);
		if (p < bound) {
		    encodedK = encoded_kmer(p);
		}
	    }
	}
    }
    if (num_hits >= min_hits) {
	process_set_of_hits(calls, otu_stats);
    }
    num_hits = 0;
}

void KmerGuts::process_aa_seq(const char *id,const char *pseq,size_t ln,
			      std::shared_ptr<std::vector<KmerCall>> calls,
			      std::function<void(sig_kmer_t &)> hit_cb,
			      std::shared_ptr<KmerOtuStats> otu_stats)
{
    strcpy(current_id,id);
    current_length_contig = ln;
    current_strand        = '+';
    current_prot_off      = 0;
    int i;
    for (i=0; (i < ln); i++)
	pIseq[i] = to_amino_acid_off(*(pseq+i));

    gather_hits(ln,'+',0,pseq,pIseq, calls, hit_cb, otu_stats);
    if (otu_stats)
	otu_stats->finalize();
}

void KmerGuts::process_seq(const char *id,const char *data,
			   std::shared_ptr<std::vector<KmerCall>> calls,
			   std::function<void(sig_kmer_t &)> hit_cb,
			   std::shared_ptr<KmerOtuStats> otu_stats)

{
    strcpy(current_id,id);
    int ln = strlen(data);
    current_length_contig = ln;
    int i;
    for (i=0; (i < 3); i++)
    {
	translate(data,i,pseq,pIseq);
	current_strand   = '+';
	current_prot_off = i;
	gather_hits(ln,'+',i,pseq,pIseq, calls, hit_cb, otu_stats);
    }
    rev_comp(data,cdata);
    for (i=0; (i < 3); i++) {
	translate(cdata,i,pseq,pIseq);

	current_strand   = '-';
	current_prot_off = i;
	gather_hits(ln,'-',i,pseq,pIseq, calls, hit_cb, otu_stats);
    }
    if (otu_stats)
	otu_stats->finalize();
}

