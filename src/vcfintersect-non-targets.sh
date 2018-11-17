#!/bin/bash

#usage: sh vcfintersect-non-targets.sh input.vcf output.vcf

vcfintersect --bed /path/to/your/bedfile/bedfile.bed --invert $1 > $2
