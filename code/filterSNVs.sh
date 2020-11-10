#! /bin/bash

# merge demultiplexed samples to one FASTQ files per barcode

# remove empty files
for f in *
do
    vcftools --vcf $f --recode --recode-INFO-all --exclude-bed /Users/gagemoreno/Desktop/2019-nCoV-tokyo.bed --out $f
done
