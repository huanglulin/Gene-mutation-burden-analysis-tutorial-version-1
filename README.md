# Gene-mutation-burden-analysis-tutorial-version-1

The rapid development of sequencing technologies has led to the launch of numerous sequencing studies for many complex traits. In addition to discovery of common variants, sequencing allows discovery of low-frequency and rare variants as well. The contribution of rare variants to disease risk is unknown for many traits, but it is reasonable to assume that rare variants influences the risk of many complex diseases. This tutorial aims to help new lab members in the Wang lab to conduct the detection of burden assays of a disease gene by using genotype data. It would tell you how to convert the genotype data and do burden assay. 
This tutorial is still being developed and we need your help to improve it: if you spot errors, or find places where clarification is needed, please let us know. If you have any suggestions to make this tutorial more useful, please also let us know. We will modify to make it easier for future readers to go through the tutorial.

0. Background
This is a background section for helping you to better understand the analysis, but skipping the background will not make the analysis not work. Thus, if you want to conduct the analysis immediately, feel free to go to the first step. Everything will work fine even without reading the background.

A. Basics of gene mutations

If you want to know more about gene mutations, please refer to https://ghr.nlm.nih.gov/primer/mutationsanddisorders/genemutation for some basic backgrounds information.

B. Basics of genotype vcf files
The VCF genotype file is very popular to store genotype data, please refer to http://support.illumina.com/content/dam/illumina-support/help/BaseSpaceHelp_v2/Content/Vault/Informatics/Sequencing_Analysis/BS/swSEQ_mBS_VCF.htm for some basic backgrounds information.

C. Basics of plink files
The Plink files is also very popular to store genotype data, please refer to http://www.cog-genomics.org/plink/1.9/ for more information.

D. Tools to be used.
Tool	version/year	Path/Link	Purpose
Plink	v.1.9	/home/lulinhuang/bin/plink	Change file format
SKAT	v 1.3.0	https://cran.r-project.org/web/packages/SKAT/	Burden assay
ANNOVAR		/home/lulinhuang/bin/  annotation

E. Training files to be used.
The training file dataset is in the following directory /share/archive/lulinhuang/geneburdentest/.


1. Variants annotation

You can make a new directory and copy the chr22_quality100.vcf file to your directory. Other files will produce during the process. 
To get the variant frequency and gene coding information for a non-annoted file, you should do annotation using ANNOVAR. For the information of ANNOVAR, please refer: http://annovar.openbioinformatics.org/en/latest/user-guide/startup/
Take a look at our chr22_quality100.vcf data:
##fileformat=VCFv4.1
##FILTER=<ID=PASS,Description="All filters passed">
##samtoolsVersion=0.1.18 (r982:295)
##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">
##INFO=<ID=DP4,Number=4,Type=Integer,Description="# high-quality ref-forward bases, ref-reverse, alt-forward and alt-reverse bases">
##INFO=<ID=MQ,Number=1,Type=Integer,Description="Root-mean-square mapping quality of covering reads">
##INFO=<ID=FQ,Number=1,Type=Float,Description="Phred probability of all samples being the same">
##INFO=<ID=AF1,Number=1,Type=Float,Description="Max-likelihood estimate of the first ALT allele frequency (assuming HWE)">
##INFO=<ID=AC1,Number=1,Type=Float,Description="Max-likelihood estimate of the first ALT allele count (no HWE assumption)">
##INFO=<ID=G3,Number=3,Type=Float,Description="ML estimate of genotype frequencies">
##INFO=<ID=HWE,Number=1,Type=Float,Description="Chi^2 based HWE test P-value based on G3">
##INFO=<ID=CLR,Number=1,Type=Integer,Description="Log ratio of genotype likelihoods with and without the constraint">
##INFO=<ID=UGT,Number=1,Type=String,Description="The most probable unconstrained genotype configuration in the trio">
##INFO=<ID=CGT,Number=1,Type=String,Description="The most probable constrained genotype configuration in the trio">
##INFO=<ID=PV4,Number=4,Type=Float,Description="P-values for strand bias, baseQ bias, mapQ bias and tail distance bias">
##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">
##INFO=<ID=PC2,Number=2,Type=Integer,Description="Phred probability of the nonRef allele frequency in group1 samples being larger (,smaller) than
in group2.">
##INFO=<ID=PCHI2,Number=1,Type=Float,Description="Posterior weighted chi^2 P-value for testing the association between group1 and group2 samples."
>
##INFO=<ID=QCHI2,Number=1,Type=Integer,Description="Phred scaled PCHI2.">
##INFO=<ID=PR,Number=1,Type=Integer,Description="# permutations yielding a smaller PCHI2.">
##INFO=<ID=VDB,Number=1,Type=Float,Description="Variant Distance Bias">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GL,Number=3,Type=Float,Description="Likelihoods for RR,RA,AA genotypes (R=ref,A=alt)">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="# high-quality bases">
##FORMAT=<ID=SP,Number=1,Type=Integer,Description="Phred-scaled strand bias P-value">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="List of Phred-scaled genotype likelihoods">
##bcftools_filterVersion=1.3.1+htslib-1.3.1
##bcftools_filterCommand=filter -e %QUAL<100 ZhenglinYang.1602EA.TruSeq.chr22.vcf
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  019032_B        019032_F        019032_M        019032_S        19070   19
071     19079   19083   19085   19086   19087   19088   19089   19090   19091   19093
.
.
.
chr22   16256352        16256352        T       C       999     PASS    DP=71;VDB=0.0298;AF1=1;AC1=2;DP4=1,0,42,25;MQ=33;FQ=-199;PV4=1,1,0.31,1 GT
:PL:GQ  0/0:0,6,28:8    0/0:0,0,0:3     0/1:16,0,8:12   0/1:5,0,26:5    0/0:0,3,12:5    0/1:7,0,20:7    0/0:0,0,0:3     0/0:0,9,29:10   0/0:0,21,6
0:22    0/0:0,15,51:16  0/0:0,5,36:7    0/1:12,0,30:11  0/0:0,84,73:80  0/1:29,0,48:28  0/0:0,3,40:5    0/1:2,0,37:4    0/0:0,0,37:4    0/1:20,0,1
0:15    0/1:25,0,25:23  0/1:13,0,27:12  0/0:0,9,29:10   0/0:0,120,55:63 0/0:0,122,57:65 0/1:9,0,41:9    0/1:37,0,18:25  0/0:0,0,44:4    0/0:0,0,46
:4      0/0:0,71,37:45  0/0:0,12,40:13  0/1:10,0,52:9   0/1:2,0,26:4    0/1:2,0,34:4    0/0:0,52,39:46  0/0:0,0,36:4    0/1:9,0,43:9    0/0:0,6,45
.
.
.


