## Generate NR clusters from BLAST database

This repo contains a Perl script used to generate non-redundant clusters of sequences from an all-vs-all BLAST database. The representatives from the resulting clusters are guaranteed to share less than X% sequence identity over Y% coverage.

This script was written a long time ago and was not designed to be run outside of the internal release pipeline in CATH. It hasn't been tested outside of CATH walls - it's here for information.

If I was doing this now, I would be tempted to look at [mmseqs](https://github.com/soedinglab/MMseqs2) or [CD-Hit](https://sites.google.com/view/cd-hit) as (much faster) alternatives.

## Install dependencies

The following seems to work in relatively modern versions of Perl (>=5.30)
 
```
cpan -I .
```

## Usage

```
$ perl script/generate_long_non_redundant_list.pl

usage: generate_long_non_redundant_list.pl [options] <S100_domain_list> <blast_results_directory>

options:

  -o|--out   <file>          Output NR file [default: nr_list_s<SEQ>_o<OV>.txt]
  --seq      Num[20-100]     Specify sequence id cutoff [default: 40.0] 
  --overlap  Num[20-100]     Specify overlap cutoff [default: 60.0] 

```
