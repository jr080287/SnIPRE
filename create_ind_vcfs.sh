
#create individual vcfs for all samples in a file

# multi-ind vcf is 1st input argument

for sample in `bcftools view -h $1 | grep "^#CHROM" | cut -f10-`; do
    bcftools view -c1 -Oz -s $sample -o $.$sample.vcf.gz $1;
done
