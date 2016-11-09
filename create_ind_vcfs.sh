#!/bin/bash

#create the individual vcfs for all samples in a file

# multi-ind vcf is 1st input argument

for sample in $(bcftools view -h  $1 | grep "^#CHROM" | cut -f10-); 
    do
        echo $sample
        bcftools view  -Oz -s $sample $1 -o $sample.vcf.gz
    done
