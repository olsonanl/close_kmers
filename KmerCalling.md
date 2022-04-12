# Using the SEED Kmer identification tools

In SEED we have two primary mechanisms for propagating function from the signature
kmer databases to input sequences.

For production use for processing protein sequences (and for making PATtyfam assignments) we
have a web service available inside on the SEED systems within the ANL firewall. The endpoint for 
the kmer function calling is `http://pear:6100/query`;