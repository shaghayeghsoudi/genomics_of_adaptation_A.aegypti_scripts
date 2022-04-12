#install.packages("tidyverse") 
library(tidyverse)
library(pophelper)
library(gridExtra)
library(mapplots)

setwd("~/Documents/malaria_vector_reserach/input_data/a.aegypti/admixture")
## samples from vcf
samp<-read.table(file = "sample_list_true_vcf.txt", header = FALSE)
samp_des<-read.table(file = "sample_list_desired.txt", header = FALSE)
rownames(samp_des)<-samp_des$V1


k1<-read.table(file = "out_k1-10/random_10K.recode.vcf.bed.1.Q", header =FALSE)
row.names(k1)<-samp[,1]

k2<-read.table(file = "out_k1-10/random_10K.recode.vcf.bed.2.Q", header =FALSE)
row.names(k2)<-samp[,1]

k3<-read.table(file = "out_k1-10/random_10K.recode.vcf.bed.3.Q", header =FALSE)
row.names(k3)<-samp[,1]

k4<-read.table(file = "out_k1-10/random_10K.recode.vcf.bed.4.Q", header =FALSE)
row.names(k4)<-samp[,1]

k5<-read.table(file = "out_k1-10/random_10K.recode.vcf.bed.5.Q", header =FALSE)
row.names(k5)<-samp[,1]

k6<-read.table(file = "out_k1-10/random_10K.recode.vcf.bed.6.Q", header =FALSE)
row.names(k6)<-samp[,1]

k7<-read.table(file = "out_k1-10/random_10K.recode.vcf.bed.7.Q", header =FALSE)
row.names(k7)<-samp[,1]

k8<-read.table(file = "out_k1-10/random_10K.recode.vcf.bed.8.Q", header =FALSE)
row.names(k8)<-samp[,1]

k9<-read.table(file = "out_k1-10/random_10K.recode.vcf.bed.9.Q", header =FALSE)
row.names(k9)<-samp[,1]

k10<-read.table(file = "out_k1-10/random_10K.recode.vcf.bed.10.Q", header =FALSE)
row.names(k10)<-samp[,1]


### order rows based on desired order
admix_K1_Q<-k1[rownames(samp_des),]
write.table(admix_K1_Q, file = "out_k1-10/out_k1-10_reordered_des/random_10K.recode.vcf.reordered.bed.1.Q", col.names= FALSE, row.names= FALSE)

admix_K2_Q<-k2[rownames(samp_des),]
write.table(admix_K2_Q, file = "out_k1-10/out_k1-10_reordered_des/random_10K.recode.vcf.reordered.bed.2.Q", col.names= FALSE, row.names= FALSE)

admix_K3_Q<-k3[rownames(samp_des),]
write.table(admix_K3_Q, file = "out_k1-10/out_k1-10_reordered_des/random_10K.recode.vcf.reordered.bed.3.Q", col.names= FALSE, row.names= FALSE)

admix_K4_Q<-k4[rownames(samp_des),]
write.table(admix_K4_Q, file = "out_k1-10/out_k1-10_reordered_des/random_10K.recode.vcf.reordered.bed.4.Q", col.names= FALSE, row.names= FALSE)

admix_K5_Q<-k5[rownames(samp_des),]
write.table(admix_K5_Q, file = "out_k1-10/out_k1-10_reordered_des/random_10K.recode.vcf.reordered.bed.5.Q", col.names= FALSE, row.names= FALSE)


admix_K6_Q<-k6[rownames(samp_des),]
write.table(admix_K6_Q, file = "out_k1-10/out_k1-10_reordered_des/random_10K.recode.vcf.reordered.bed.6.Q", col.names= FALSE, row.names= FALSE)


admix_K7_Q<-k7[rownames(samp_des),]
write.table(admix_K7_Q, file = "out_k1-10/out_k1-10_reordered_des/random_10K.recode.vcf.reordered.bed.7.Q", col.names= FALSE, row.names= FALSE)


admix_K8_Q<-k8[rownames(samp_des),]
write.table(admix_K8_Q, file = "out_k1-10/out_k1-10_reordered_des/random_10K.recode.vcf.reordered.bed.8.Q", col.names= FALSE, row.names= FALSE)


admix_K9_Q<-k9[rownames(samp_des),]
write.table(admix_K9_Q, file = "out_k1-10/out_k1-10_reordered_des/random_10K.recode.vcf.reordered.bed.9.Q", col.names= FALSE, row.names= FALSE)


admix_K10_Q<-k10[rownames(samp_des),]
write.table(admix_K10_Q, file = "out_k1-10/out_k1-10_reordered_des/random_10K.recode.vcf.reordered.bed.10.Q", col.names= FALSE, row.names= FALSE)


####### plot the results
## admixture format 

setwd("~/Documents/malaria_vector_reserach/input_data/a.aegypti/admixture/out_k1-10")
afiles  <- list.files("out_k1-10_reordered_des", pattern="*Q", full.names = TRUE)
alist <- readQ(files=afiles)


attributes(alist[[1]])


tr1 <- tabulateQ(qlist=alist)


tabulateQ(alist)
tabulateQ(alist, writetable=TRUE)
sr1 <- summariseQ(tr1)


p1 <- plotQ(alist[c(1,3,4,5,6,7,8,9,2)],imgoutput="join",returnplot=T,exportplot=F,quiet=T,basesize=11)
grid.arrange(p1$plot[[1]])



p1<-plotQ(alist[c(1,3,4,5,6)], imgoutput="join",splab=spnames[c(1,3,4,5,6)] ,splabcol="blue",splabsize=9,returnplot=T,exportplot=F,quiet=T,basesize=1,showyaxis=T,grplab=onelabset1,,grplabsize=3,linesize=0.7,pointsize=5)
grid.arrange(p1$plot[[1]])


pops <- read.delim(file ="~/Documents/malaria_vector_reserach/input_data/a.aegypti/admixture/metadata_orderd.txt",header=FALSE,stringsAsFactors=F)
colnames(pops)<-c("ind","region")

fn1 <- function(x) attr(x,"k")
spnames <- paste0("K=",sapply(alist,fn1))
onelabset1 <- pops[,2,drop=FALSE]
head(onelabset1)
p1<-plotQ(qlist=alist, imgoutput="join",splab=spnames[1:10],splabcol="blue",splabsize=7,returnplot=T,exportplot=T,quiet=T,basesize=20,showyaxis=T,indlabsize=5,grplab=onelabset1,grplabsize=3,linesize=0.7,pointsize=5)
grid.arrange(p1$plot[[1]])





#!/usr/bin/Rscript

snps<-read.delim('results.012',header=F,row.names=1)
pos<-read.delim('results.012.pos',header=F)
indv<-read.delim('results.012.indv',header=F)
colnames(snps)<-paste(pos[,1],pos[,2],sep=':')
rownames(snps)<-indv[,1]
snps<-as.matrix(snps)
snps.convert<-t(snps)

write.table(snps.convert, file= "snps.convert_for_lfmm_inds", col.names = FALSE, row.names = TRUE)

