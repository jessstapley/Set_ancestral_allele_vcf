# Setting the ancestral allele in a vcf file for haploid and diploid data
This document describes how to annotate a vcf file to indicate which allele is the ancestral allele. 

I still need to finish the scripts that it references for Step 1-4 - download sequence reads from the NCBI Short Read Arcives (SRA), how to map reads to the reference genome and to call variants using the GATK pipeline.

You can begin from Step 6 if you have a vcf file that includes all your samples and the outgroup. The appraoch closely follows instructions posted here 
http://wasabiapp.org/vbox/data/ngg2016/21/Day2Session1.Genomicalignmentancestralalleles.html

Here I provide instructions on how to do this for both Hapoid and Diploid data. You need a basic working knowledge of the command line and experience with vcftools and bcftools would be beneficial.

# What you need
A draft genome sequence for your focal species.
Sequence reads (e.g. from SRA) from another closely related species - an outgroup.

If you have sequence reads available from another species that is closely related to your focal species then you can map these reads to the reference genome of your focal species, call variants and determine the ancestral state for the variants.

NB: How genetically close your outgroup is and how well your outgroup maps to your reference will largely determine the number of sites where the ancestral allele can be verified.

# Software Requirements
SRA manipulation:
SRAtools
fastqc
trimmomatic

For mapping and variant calling: 
bwa mem
GATK
vcftools

for setting ancestral allele from vcf file:
samtools

# Step by step guide

Step 1. Extract the sequence reads from SRA, check the quality, remove adaptors and trim (see here)

Step 2. Map with bwa mem (see here)

Step 3. Call variants with GATK (see here) - all indviduals combined
	
Step 4. Filter variants, remove indels and non-biallelic SNPs using vcftools
e.g biallelic SNPs only
vcftools --vcf file.vcf  --max-alleles 2 --recode --recode-INFO-all --out file_BI
e.g SNPs only
vcftools --vcf file_BI.recode.vcf  ---remove-indels --recode --recode-INFO-all --out file_BI_SNPS


Step 5. Create a vcf file for the outgroup individual only
vcftools --vcf file_BI_SNPS.recode.vcf --indv outgroupIDname --recode --recode-INFO-all --out outgroupIDname

Step 6. Create a table of SNP positions and alleles using bcftools 
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' outgroupIDname.recode.vcf  > file.tab

Step 7. Create a table file with ancestral allele information.
Here we select the reference allele (REF) or the variant allele (ALT) and put it in the 5th column. If the site is not covered by our outgroup reads, i.e. missing '.', then we record it as missing.

for Haploid    
awk '{OFS="\t";if($5=="0"){print $1,$2,$3,$4,$3} \
	if($5=="1"){print $1,$2,$3,$4,$4} \
	if($5=="."){print $1,$2,$3,$4,$5}}' file.tab > file_aa.tab


for Diploid
awk '{OFS="\t";if($5=="0/0"){print $1,$2,$3,$4,$3} \
	if($5=="0/1"){print $1,$2,$3,$4,$4} \
	if($5=="./."){print $1,$2,$3,$4,$5}}' file.tab > file_aa.tab
                  
   
Step 8. Compress and index the table file and the original vcf file
bgzip file_aa.tab
tabix -s1 -b2 -e2 file_aa.tab.gz
bgzip file_BI_SNPS.recode.vcf

Step 9. Create an INFO file line for the new vcf file
echo '##INFO=<ID=AA,Number=1,Type=Character,Description="Ancestral allele">' > hdr.txt


Step 10. Using bcftools to annotate the vcf file with the ancestral allele information 
bcftools annotate -a file_aa.tab.gz \
 -c CHROM,POS,REF,ALT,INFO/AA -h hdr.txt -Oz \
 -o newfile_aa.vcf.gz file_BI_SNPS.recode.vcf.gz


Step 11. Check that it has worked. There should be an info field AA
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AA\n' newfile_aa.vcf.gz | less
               
Step 12. Count how many sites have ancestral allele information
bcftools view -e 'INFO/AA=="."' newfile_aa.vcf.gz -H | wc -l

