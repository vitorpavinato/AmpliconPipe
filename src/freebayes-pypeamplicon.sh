#!/bin/bash 

#usage: sh freebayes-pypeamplicon.sh input.bam output.vcf populations-file

freebayes -b $1 -f /path/to/your/reference/amplicons.fa -t /path/to/your/bedfile/amplicons_regions.bed -v $2 --report-all-haplotype-alleles --strict-vcf --ploidy 2 --use-best-n-alleles 6 -E 9 -4 --populations $3 --genotype-qualities

