#!/usr/bin/env Rscript
## PCA analysis bu SNp-relate with entire data set_H.annuus_2016_GWS ##

#Load the pakcages gdsfmt and SNPRelate and tidyverse(for plotting) into your R session #
#library(gdsfmt)
library(SNPRelate)

## Set your working directory to where the files are located
setwd("/data/home/shaghayegh/vgl/test_ag/maf_0.05/ld_snp_relate")


######################################################################################################
## prepare input file directly from vcf file

## give path to the vcf file (my working directory is the folder that the vcf file is located)
## I have filtered vcf file for 315 good individuals and filtered DP15 SNPs

vcf_filename<- c("/data/home/shaghayegh/vgl/test_ag/maf_0.05/YL-Agam-GF2_pflit.soft.mask.good.inds.recode.biall.indel.removed.Q30.DP6.maxmissing0.8.maf0.05.recode.vcf")

## reformat
snpgdsVCF2GDS(vcf_filename, "YL-Agam-GF2_pflit.soft.mask.good.inds.recode.biall.indel.removed.Q30.DP6.maxmissing0.8.maf0.05.gds", method="biallelic.only")

## get summary
snpgdsSummary("YL-Agam-GF2_pflit.soft.mask.good.inds.recode.biall.indel.removed.Q30.DP6.maxmissing0.8.maf0.05.gds")

## to open gds file
(genofile <- snpgdsOpen("YL-Agam-GF2_pflit.soft.mask.good.inds.recode.biall.indel.removed.Q30.DP6.maxmissing0.8.maf0.05.gds"))


#sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
## make all info created by snprelate
snp.chrom <- read.gdsn(index.gdsn(genofile, "snp.chromosome"))
snp.chrom.dat<-data.frame(snp.chrom)

snp.pos <- read.gdsn(index.gdsn(genofile, "snp.position"))
snp.pos.dat<-data.frame(snp.pos)

snp.id <- read.gdsn(index.gdsn(genofile, "snp.id"))
snp.id.dat<-data.frame(snp.id)

dat.snp.all.info<-cbind(snp.chrom.dat,snp.pos.dat,snp.id.dat)
write.table(dat.snp.all.info,file = "YL-Agam-GF2_pflit.soft.mask.good.inds.recode.biall.indel.removed.Q30.DP6.maxmissing0.8.maf0.05.SNPrelate.info",col.names = TRUE)






## data analysis ##
## Linkage Disequilibrium (LD) Based SNP Pruning
snpset_pruned <- snpgdsLDpruning(genofile,ld.threshold=0.1, slide.max.bp = 2000, maf = 0.05, autosome.only=F)

#Excluding 12,457 SNPs (monomorphic: TRUE, < MAF: NaN, or > missing rate: NaN)
#Working space: 315 samples, 11,071,304 SNPs
#    using 1 (CPU) core
#	Sliding window: 50000 basepairs, Inf SNPs
#	|LD| threshold: 0.1
#Chromosome HanXRQChr01: 5.89%, 35392/601001
#Chromosome HanXRQChr02: 6.28%, 34153/543849
#Chromosome HanXRQChr03: 5.60%, 35628/635995
#Chromosome HanXRQChr04: 5.92%, 38023/642748
#Chromosome HanXRQChr05: 6.11%, 44264/724770
#Chromosome HanXRQChr06: 5.55%, 23179/417530
#Chromosome HanXRQChr07: 5.66%, 22675/400390
#Chromosome HanXRQChr08: 5.71%, 35029/613066
# ...

names(snpset_pruned)

##  get a list of prunned SNPs
snpset.id <- unlist(snpset_pruned)


snpset.id.g <- data.frame(unlist(snpset_pruned))
write.table(snpset.id.g, file = "YL-Agam-GF2_pflit.soft.mask.good.inds.recode.biall.indel.removed.Q30.DP6.maxmissing0.8.maf0.05.SNPrelate.slide2000.for_str",col.names = FALSE, row.names = FALSE)

overlap<- dat.snp.all.info[(dat.snp.all.info$snp.id %in% snpset.id.g$V1),]
write.table(overlap, file = "YL-Agam-GF2_pflit.soft.mask.good.inds.recode.biall.indel.removed.Q30.DP6.maxmissing0.8.maf0.05.SNPrelate.slide2000.for_str_LD_prunned.txt", row.names = FALSE, col.names = True)


