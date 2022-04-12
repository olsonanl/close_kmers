# Using the SEED Kmer identification tools

In SEED we have two primary mechanisms for propagating function from the signature
kmer databases to input sequences.

## Kmer calling web service

For production use for processing protein sequences (and for making PATtyfam assignments) we
have a web service available inside on the SEED systems within the ANL firewall. The endpoint for 
the kmer function calling is `http://pear:6100/query`. The protein sequence is passed as the input 
to the HTTP POST request. For example, if the file query.fa holds one or more protein sequences in FASTA format
the curl command may be used to process the data:

```
$ curl --data-binary @query.fa http://pear:6100/query
PROTEIN-ID	     id-83333.1.peg.1	21
CALL		     0			19	13	30895	Thr operon leader peptide	22.729
OTU-COUNTS	     id-83333.1.peg.1[21]	1-91	1-200	1-251	   1-460  1-475
PROTEIN-ID	     id-83333.1.peg.2		820
CALL		     0				818	791	3449	Aspartokinase (EC 2.7.2.4) / Homoserine dehydrogenase (EC 1.1.1.3)	2432.37
OTU-COUNTS	     id-83333.1.peg.2[820]	8-310	7-398	6-163	6-244	      5-11
PROTEIN-ID	     id-83333.1.peg.3		310
CALL		     0				307	293	15456	Homoserine kinase (EC 2.7.1.39)	804.771
OTU-COUNTS	     id-83333.1.peg.3[310]	4-907	4-208	3-189	3-695	   3-378
PROTEIN-ID	     id-83333.1.peg.4		428
CALL		     0				426	412	30907	Threonine synthase (EC 4.2.3.1)	1232.38
OTU-COUNTS	     id-83333.1.peg.4[428]	8-310	5-438	5-623	4-3	  4-12
```

The fields in the output are as follows.

 * PROTEIN-ID: identifier (from fasta file), sequence length

 * CALL: start of called region (in the protein sequence), end of called region, count of kmer hits

 * OTU-COUNTS: sequence id, sequence length (in brackets), and a list of otu signatures of the format OTU-Count

## `kmer_search` command line tool