#!/bin/bash
#
#SBATCH --time=10:00:00
#SBATCH --qos=normal
#SBATCH --ntasks=1
#SBATCH --mem=8000

#usage: sh vcfintersect-targets.sh input.vcf output.vcf

vcfintersect --bed /home/correapavinato1/Working_dir/Data/Fluidigm/Belgica/bedfile/Belgica_target_SNP_region.bed $1 > $2
