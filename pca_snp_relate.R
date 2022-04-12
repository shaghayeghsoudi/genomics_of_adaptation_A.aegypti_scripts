#!/usr/bin/env Rscript
rm(list=ls())
## PCA analysis bu SNp-relate with entire data set_H.annuus_2016_GWS ##

#Load the pakcages gdsfmt and SNPRelate and tidyverse(for plotting) into your R session #
#library(gdsfmt)
library(SNPRelate)
#display.brewer.all()
library(RColorBrewer)
library(ggplot2)


## Set your working directory to where the files are located
setwd("/share/lanzarolab/users/shaghayegh/data/malli_vs_benin/data_March2020/pca/2L_2-3MB_48inds_phylo")


######################################################################################################
## prepare input file directly from vcf file

## give path to the vcf file (my working directory is the folder that the vcf file is located)
## I have filtered vcf file for 315 good individuals and filtered DP15 SNPs

vcf_filename<- c("2L_kdr_hom.recode.vcf")

## reformat
snpgdsVCF2GDS(vcf_filename, "2L_2-3M.gds", method="biallelic.only")

## get summary
snpgdsSummary("2L_2-3M.gds")

## to open gds file
(genofile <- snpgdsOpen("2L_2-3M.gds"))


#sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
## make all info created by snprelate
snp.chrom <- read.gdsn(index.gdsn(genofile, "snp.chromosome"))
snp.chrom.dat<-data.frame(snp.chrom)

snp.pos <- read.gdsn(index.gdsn(genofile, "snp.position"))
snp.pos.dat<-data.frame(snp.pos)

snp.id <- read.gdsn(index.gdsn(genofile, "snp.id"))
snp.id.dat<-data.frame(snp.id)

dat.snp.all.info<-cbind(snp.chrom.dat,snp.pos.dat,snp.id.dat)
write.table(dat.snp.all.info,file = "2L.2-3MB.pos.48.inds.phylo.gds.info.txt",col.names = TRUE, quote= FALSE)


## data analysis ##
## Linkage Disequilibrium (LD) Based SNP Pruning
snpset_pruned <- snpgdsLDpruning(genofile,ld.threshold=0.2, slide.max.bp = 1000, autosome.only=F)

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
colnames(snpset.id.g)<-c("id")
write.table(snpset.id.g, file = "2L.pos_1_3MB.gds.SNPrelate.slide1000.ld.threshold.0.2.for_str.txt",col.names = FALSE, row.names = FALSE, quote = FALSE)

overlap<- dat.snp.all.info[(dat.snp.all.info$snp.id %in% snpset.id.g$id),]
write.table(overlap, file = "2L.2-3MB.gds.slide1000.ld.threshold.0.2.for_str_LD_prunned.txt", row.names = FALSE, col.names = TRUE,quote = FALSE)

######################
## do PCA analysis ##
######################
pca <- snpgdsPCA(genofile, num.thread = 2, snp.id = snp.id,autosome.only = F)
pc.percent <- pca$varprop*100
pc.round<-head(round(pc.percent, 2))
write.table(pc.round,file = "2L.2-3MB.pos.48.inds.phylo_variance_proportion_pca.txt",quote = FALSE)

################################################
#### if there is no population information
################################################

## when there is no population information

tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2], 
                  # the second eigenvector
                  stringsAsFactors = FALSE)
head(tab)

#sample.id         EV1         EV2
#1   HT185_MK8-2.sort.bam -0.03885285 0.047723358
#2  HT48_SD1A-37.sort.bam  0.04219809 0.031562977
#3  HT53_IA1A-10.sort.bam -0.02669614 0.003319497
#4   HT260_MK1-4.sort.bam -0.06099702 0.008044707
#5 HT170_SD2W-14.sort.bam  0.01353301 0.010424208
#6   Q093_MK7-13.sort.bam  0.06154042 0.155766071


###### get population info ######
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

# Get population information
pop_code <- scan("pop.txt", what=character())
#   if it is stored in a text file "pop.txt"

head(cbind(sample.id, pop_code))
#sample.id                pop_code
#[1,] "HT185_MK8-2.sort.bam"   "HA_27"
#[2,] "HT48_SD1A-37.sort.bam"  "HA_14"
#[3,] "HT53_IA1A-10.sort.bam"  "HA_1"
#[4,] "HT260_MK1-4.sort.bam"   "HA_20"
#[5,] "HT170_SD2W-14.sort.bam" "HA_17"
#[6,] "Q093_MK7-13.sort.bam"   "HA_26"


tab <- data.frame(sample.id = pca$sample.id,
                  pop = factor(pop_code)[match(pca$sample.id, sample.id)],
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  EV3 = pca$eigenvect[,3],
                  EV4 = pca$eigenvect[,4],
                  EV5 = pca$eigenvect[,5],
                  EV6 = pca$eigenvect[,6],
                  stringsAsFactors = FALSE)
head(tab)

#sample.id   pop         EV1         EV2
#1   HT185_MK8-2.sort.bam HA_27 -0.03885285 0.047723358
#2  HT48_SD1A-37.sort.bam HA_14  0.04219809 0.031562977
#3  HT53_IA1A-10.sort.bam  HA_1 -0.02669614 0.003319497
#4   HT260_MK1-4.sort.bam HA_20 -0.06099702 0.008044707
#5 HT170_SD2W-14.sort.bam HA_17  0.01353301 0.010424208
#6   Q093_MK7-13.sort.bam HA_26  0.06154042 0.155766071

write.table(tab, file = "pca_variance_proportionpconcat_2L.2-3MB.pos.48.inds.phylo.txt")



###########################################################################
## plot PC plot
### we have 71 populations, to adjust the color
## tutorial : 
## https://www.r-graph-gallery.com/colors/
#library(ggplot2)
#library(RColorBrewer)
setwd("~/Documents/malaria_vector_reserach/input_data/introgression_mali_benin/data_March_2020/2-3MB_kdr_hom/pca")
pca<-read.table(file = "pca_variance_proportionpconcat_2L.2-3MB.pos.48.inds.phylo.txt", header = TRUE)

pdf("2L_2-3MB_pos_all_inds_2020_vcf_PC1-PC2.pdf", height= 6, width = 6)
print(ggplot(pca,aes(x=EV1,y=EV2,color=pop))+geom_point(size=5,alpha = .5)+scale_size_area()+scale_fill_brewer(palette="Set3"))
dev.off()


pdf("2L_2-3MB_pos_all_inds_2020_vcf_PC3-PC4.pdf",height= 6, width = 6)
print(ggplot(pca,aes(x=EV3,y=EV4,color=pop))+geom_point(size=5,alpha = .5)+scale_size_area()+scale_fill_brewer(palette="Set3"))
dev.off()

pdf("2L_2-3MB_pos_all_inds_2020_vcf_PC5-PC6.pdf",height= 6, width = 6)
print(ggplot(pca,aes(x=EV5,y=EV6,color=pop))+geom_point(size=5,alpha = .5)+scale_size_area()+scale_fill_brewer(palette="Set3"))
dev.off()
