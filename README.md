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
##INFO=<ID=PC2,Number=2,Type=Integer,Description="Phred probability of the nonRef allele frequency in group1 samples being 
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


In this case, you can use the following command to convert the vcf file to vcf4 for ANNOVAR annotation:

[lulinhuang@compute-0-54 geneburdentest]$ convert2annovar.pl -format vcf4 chr22_quality100.vcf -outfile vcf4chr22_quality100

Then, using ANNOVAR annotation to get the gene, variant frequency and mutation CADD score for further analysis.

[lulinhuang@compute-0-54 geneburdentest]$ table_annovar.pl vcf4chr22_quality100 -thread 6 -buildver hg19 /home/lulinhuang/bin/annovar/humandb -out vcf4chr22_quality100 -remove -protocol refGene,cytoBand,1000g2015aug_all,cadd13 -operation g,r,f,f -nastring .

These commands will give the annotation file:

[lulinhuang@compute-0-54 geneburdentest]$ more vcf4chr22_quality100_SNPs.hg19_multianno.txt

Chr     Start   End     Ref     Alt     Func.refGene    Gene.refGene    GeneDetail.refGene      ExonicFunc.refGene      AAChange.refGene        cytoBand  1000g2015aug_all        CADD13_RawScore CADD13_PHRED
chr22   16269934        16269934        A       G       exonic  POTEH   .       nonsynonymous SNV       POTEH:NM_001136213:exon7:c.T1247C:p.L416S22q11.1  .       2.690721        20.8
chr22   16287789        16287789        C       G       exonic  POTEH   .       nonsynonymous SNV       POTEH:NM_001136213:exon1:c.G97C:p.A33P  22q11.1   0.0948482       -0.763275       0.052

2. Make the input files

2.1 Make binary plink files

Using the following commands to convert the vcf genotype file to get the binary plink data files which will be used as part of the input files for SKAT package.

[lulinhuang@compute-0-54 geneburdentest]$ plink --vcf chr22_quality100.vcf --make-bed --double-id --out chr22_quality100

This command will get these three files:

-rw-r--r--  1 lulinhuang wanglab   5712248 Jul 24 19:52 chr22_quality100.bed
-rw-r--r--  1 lulinhuang wanglab    431279 Jul 24 19:52 chr22_quality100.bim
-rw-r--r--  1 lulinhuang wanglab     33977 Jul 24 19:52 chr22_quality100.fam

You can also do association analysis using these files to get the P values to do QQ plot to see the data quality.

The association analysis command:

[lulinhuang@compute-0-54 geneburdentest]$ plink --bfile chr22_quality100 --assoc --out chr22_quality100

[lulinhuang@compute-0-54 geneburdentest]$ more chr22_quality100.assoc

 CHR           SNP         BP   A1      F_A      F_U   A2        CHISQ            P           OR
  22   22:16256352   16256352    C   0.2681   0.2796    T      0.08198       0.7746       0.9439
  22   22:16256430   16256430    G   0.1812   0.1723    A      0.06807       0.7942        1.063
  22   22:16256512   16256512    C   0.2609   0.2463    T       0.1426       0.7057         1.08
  22   22:16269934   16269934    G   0.4203   0.4297    A      0.04559       0.8309       0.9621

Then, you can use qqman package to do QQ plot using the association data.

install.packages("qqman")

>library(qqman)

qq(chr22_quality100.assoc$P, main = "Q-Q plot of GWAS p-values")

Then, you will got the Q-Q plot like this:

2.2 Make File.SetID file

Extract the variants with the frequency minor allele frequency (MAF) in 1000g2015aug_all database less than 0.01 or MAF less than 0.001 from the annotation file vcf4chr22_quality100_SNPs.hg19_multianno.txt, because we want to analysis rare variants.

The gene setID file is like this:

[lulinhuang@compute-0-54 geneburdentest]$ more chr22_quality100_0.01_missense.SetID

POTEH   22:16277873
POTEH   22:16277880
POTEH   22:16279292
POTEH   22:16287354
POTEH   22:16287495
POTEH   22:16287753
POTEH   22:16287794
OR11H1  22:16449210
OR11H1  22:16449320
CCT8L2  22:17071802
CCT8L2  22:17072285
CCT8L2  22:17072384
CCT8L2  22:17072486
CCT8L2  22:17072504

[lulinhuang@compute-0-54 geneburdentest]$ more chr22_quality100_0.001_missense.SetID

POTEH   22:16279292
POTEH   22:16287495
POTEH   22:16287753
POTEH   22:16287794
CCT8L2  22:17071802
CCT8L2  22:17072285
CCT8L2  22:17072486
CCT8L2  22:17072504
CCT8L2  22:17072743
CCT8L2  22:17072806
CCT8L2  22:17073028
CCT8L2  22:17073131
CCT8L2  22:17073133

2.3 Make Weight file using the CADD score of the ANNOVAR annotation file vcf4chr22_quality100_SNPs.hg19_multianno.txt (CADD13_PHRED score, the larger value indicate more functional damage for a variant). For more information about the CADD score, please infer http://cadd.gs.washington.edu/info. 

The Weight file is like this:

[lulinhuang@compute-0-54 geneburdentest]$ more  missense_0.01_chr22_quality100weights

22:16277873     12.54
22:16277880     19.3
22:16279292     23.4
22:16287354     17.25
22:16287495     0.008
22:16287753     0.001
22:16287794     18.9
22:16449210     35
22:16449320     0.192
22:17071802     0.001
22:17072285     23.1
22:17072384     0.053
22:17072486     23

[lulinhuang@compute-0-54 geneburdentest]$ more missense_0.001_chr22_quality100weights

22:16279292     23.4
22:16287495     0.008
22:16287753     0.001
22:16287794     18.9
22:17071802     0.001
22:17072285     23.1
22:17072486     23
22:17072504     0.68
22:17072743     13.38
22:17072806     1.615
22:17073028     5.251
22:17073131     25.9
22:17073133     24
22:17073145     18.72
22:17073406     22.4
22:17288788     19.19
22:17444640     15.02
22:17444685     23.8
22:17445711     6.354

3. Running SKAT 

Install SKAT R package:
[lulinhuang@compute-0-54 geneburdentest]$R
> install.packages("SKAT_1.3.0.tar.gz", repos = NULL, type="source")

For running the missense maf 0.01 dataset, using the following R scripts:

library(SKAT)
setwd("/share/archive/lulinhuang/geneburdentest")
Project.BED="chr22_quality100.bed"
Project.BIM="chr22_quality100.bim"
Project.FAM="chr22_quality100.fam"
Project.SetID="chr22_quality100_0.01_missense.SetID"
Project.SSD="chr22.SSD"
Project.Info="chr22.SSD.info"
Generate_SSD_SetID(Project.BED,Project.BIM,Project.FAM,Project.SetID,Project.SSD,Project.Info)
Project.FAM="chr22_quality100.fam"
FAM<-Read_Plink_FAM(Project.FAM, Is.binary=FALSE)
y<-FAM$Phenotype
Project.SSD="chr22.SSD"
Project.Info="chr22.SSD.info"
SSD.INFO<-Open_SSD(Project.SSD,Project.Info)
weights="missense_0.01_chr22_quality100weights"
Read_SNP_WeightFileweights="weights"
mw = Read_SNP_WeightFile(weights)
obj<-SKAT_Null_Model(y ~ 1, out_type="C")
missense_out0.01chr22=SKAT.SSD.All(SSD.INFO,obj,method="davies",obj.SNPWeight = mw)
write.csv(as.data.frame(missense_out0.01chr22$results),file=paste("./missense_out0.01chr22_davies",".csv", sep=""))

The results are like this:

[lulinhuang@compute-0-54 geneburdentest]$ more missense_out0.01chr22_davies.csv

"","SetID","P.value","N.Marker.All","N.Marker.Test"
"1","3-Sep",0.571489184186594,4,2
"2","5-Sep",1,3,3
"3","A4GALT",0.97498862704667,4,4
"4","ACO2",0.0754998779690406,5,3
"5","ACR",0.68719124428967,3,3
"6","ADM2",0.611640451068496,1,1
"7","ADORA2A",0.849058913694201,5,4

For running the missense maf 0.001 dataset, using the following R scripts:

library(SKAT)
setwd("/share/archive/lulinhuang/geneburdentest")
Project.BED="chr22_quality100.bed"
Project.BIM="chr22_quality100.bim"
Project.FAM="chr22_quality100.fam"
Project.SetID="chr22_quality100_0.001_missense.SetID"
Project.SSD="chr22.SSD"
Project.Info="chr22.SSD.info"
Generate_SSD_SetID(Project.BED,Project.BIM,Project.FAM,Project.SetID,Project.SSD,Project.Info)
Project.FAM="chr22_quality100.fam"
FAM<-Read_Plink_FAM(Project.FAM, Is.binary=FALSE)
y<-FAM$Phenotype
Project.SSD="chr22.SSD"
Project.Info="chr22.SSD.info"
SSD.INFO<-Open_SSD(Project.SSD,Project.Info)
weights="missense_0.001_chr22_quality100weights"
Read_SNP_WeightFileweights="weights"
mw = Read_SNP_WeightFile(weights)
obj<-SKAT_Null_Model(y ~ 1, out_type="C")
missense_out0.001chr22=SKAT.SSD.All(SSD.INFO,obj,method="davies",obj.SNPWeight = mw)
write.csv(as.data.frame(missense_out0.001chr22$results),file=paste("./missense_out0.001chr22_davies",".csv", sep=""))

The results are like this:

[lulinhuang@compute-0-54 geneburdentest]$ more missense_out0.001chr22_davies.csv

"","SetID","P.value","N.Marker.All","N.Marker.Test"
"1","3-Sep",0.769654727591407,3,1
"2","5-Sep",1,3,3
"3","A4GALT",0.97498862704667,4,4
"4","ACO2",0.769654727591407,3,1
"5","ACR",0.68719124428967,3,3
"6","ADORA2A",0.849058918816583,4,3
"7","ADRBK2",1,3,2
"8","ADSL",0.0450475080174264,3,2
"9","AIFM3",0.061376924737658,11,7
"10","ALG12",0.769713635200432,5,2
"11","ANKRD54",1,2,2
"12","AP1B1",0.138144373372949,4,2
"13","APOBEC3A_APOBEC3A_B",0.0305064455513837,3,2
"14","APOBEC3B",0.769654724897054,3,2
"15","APOBEC3C",0.918974989823189,8,6
"16","APOBEC3D",0.0331526980626367,5,3
"17","APOBEC3F",0.969232906360845,5,3
"18","APOBEC3G",0.967268302448477,7,5
"19","APOBEC3H",0.769654720596527,2,2
"20","APOL1",0.942743752448574,8,7
"21","APOL2",0.781510937878202,5,4






