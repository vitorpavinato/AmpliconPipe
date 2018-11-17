#!/bin/bash

#usage: sh vcfintersect-targets.sh input.vcf output.vcf

vcfintersect --bed /path/to/your/bedfile/bedfile.bed $1 > $2
