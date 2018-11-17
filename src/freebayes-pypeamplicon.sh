#!/bin/bash 
#
#SBATCH --time=10:00:00
#SBATCH --qos=normal
#SBATCH --ntasks=1
#SBATCH --mem=8000

#usage: sh freebayes-ampliconseq-belgica.sh input.bam output.vcf populations-file

freebayes -b $1 -f /home/correapavinato1/Working_dir/Data/Fluidigm/Belgica/reference/Belgica_fluidigm_amplicons.fa -t /home/correapavinato1/Working_dir/Data/Fluidigm/Belgica/bedfile/Belgica_amplicons_regions.bed -v $2 --report-all-haplotype-alleles --strict-vcf --ploidy 2 --use-best-n-alleles 6 -E 9 -4 --populations $3 --genotype-qualities

