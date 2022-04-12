#!/usr/bin/env Rscript
library("qvalue")
setwd("/share/lanzarolab/users/shaghayegh/data/A.aegypti/good_ind_Q30.DP5.caling_1_maf0.03/downstream_analysis/lfmm")
all_good_chrom_lfmm<-read.table(file = "var_out_A.aegypti_Q30.DP5.caling_1_maf0.03_lfmm_IS_STD_12_pops_with_genes_attached.txt", header = FALSE)
colnames(all_good_chrom_lfmm)<-c("snp_id","chr","pos","latitude_e","longitude_e","MAT_e","MWMT_e","MCMT_e","TD_e","MAP_e","MSP_e","AHM_e","SHM_e","DD_0_e","DD5_e","DD_18_e","DD18_e","NFFD_e","bFFP_e","eFFP_e","FFP_e","PAS_e","EMT_e","EXT_e","Eref_e","CMD_e","MAR_e","RH_e","genes")


test_type_env<-array ("0",ncol (all_good_chrom_lfmm))
env_type<-grep ("_e$",colnames(all_good_chrom_lfmm))
test_type_env[min(env_type):(min(env_type)+24)]<-"envir"

all_good_chrom_lfmm$id<- paste (all_good_chrom_lfmm[,2],all_good_chrom_lfmm[,3], sep = "__")

out_res <- NULL
q = 0.01

test_names_env <- colnames (all_good_chrom_lfmm)[test_type_env == "envir"]
the_i <- which (test_type_env == "envir") 



for (i in 1:length (the_i)){
  L <- as.numeric(length(all_good_chrom_lfmm[,the_i[i]]))
  w = which(sort(all_good_chrom_lfmm[,the_i[i]]) < q * (1:L) / L)
  candidates = all_good_chrom_lfmm[order(all_good_chrom_lfmm[,the_i[i]])[w],c(1,2,3,the_i[i],29)]
  candidates$var<-rep(test_names_env[i],nrow(candidates))
  write.table(candidates, file = paste("top_candidate_genes_FDR0.01/lfmm_A.aegypti_12_pops_list_candidate_loci_FDR_0.01_with_attached_genes_vars_",test_names_env[i],".txt", sep = ""),col.names = TRUE, row.names= FALSE,quote= FALSE)
}


#######################



#!/usr/bin/env Rscript
library(qqman)

setwd("/share/lanzarolab/users/shaghayegh/data/A.aegypti/good_ind_Q30.DP5.caling_1_maf0.03/downstream_analysis/lfmm")
#setwd("~/Documents/malaria_vector_reserach/input_data/a.aegypti/lfmm")


all_good_chrom_lfmm<-read.table(file="var_out_A.aegypti_Q30.DP5.caling_1_maf0.03_lfmm_IS_STD_12_pops_with_genes_attached.txt", header = FALSE)
#all_good_chrom_lfmm<-read.table(file="var_out_lfmm_test.txt", header = FALSE)

colnames(all_good_chrom_lfmm)<-c("snp_id","chr","pos","latitude_e","longitude_e","MAT_e","MWMT_e","MCMT_e","TD_e","MAP_e","MSP_e","AHM_e","SHM_e","DD_0_e","DD5_e","DD_18_e","DD18_e","NFFD_e","bFFP_e","eFFP_e","FFP_e","PAS_e","EMT_e","EXT_e","Eref_e","CMD_e","MAR_e","RH_e","genes")

#all_lfmm_genes <- all_good_chrom_lfmm[grep("gene",all_good_chrom_lfmm$genes),]
test_type_env<-array ("0",ncol (all_good_chrom_lfmm))
env_type<-grep ("_e$",colnames(all_good_chrom_lfmm))
test_type_env[min(env_type):(min(env_type)+24)]<-"envir"

all_good_chrom_lfmm$id<- paste (all_good_chrom_lfmm[,2],all_good_chrom_lfmm[,3], sep = "__")

out_res <- NULL
q = 0.01

test_names_env <- colnames (all_good_chrom_lfmm)[test_type_env == "envir"]
the_i <- which (test_type_env == "envir") 

genes<-read.table(file = "lfmm_A.aegypti_12_pops_list_candidate_loci_FDR_0.01_with_attached_genes_vars_ALL_Vars.txt", header = TRUE)
colnames(genes)<-c("snp_id","CHR","pos","BF_is","genes","var")
chromosome<-unique(all_good_chrom_lfmm$CHR)


for (i in 1:length (the_i)){
  
  genes_g<-data.frame("gene"=genes[genes$var==test_names_env[i],5])
  genes_g<-na.omit(genes_g)
  #genes_zoom<-data.frame("gene"=genes_g[genes_g$chrom==j,2])
  all_good_chrom_a<-all_good_chrom_lfmm[,c(1,2,3,the_i[i],29)]
  
  
  #L <- as.numeric(length(all_good_chrom_a[,the_i[i]]))
  #w = which(sort(all_good_chrom_a[,the_i[i]]) < q * (1:L) / L)
  #candidates = all_good_chrom_lfmm[order(all_good_chrom_lfmm[,the_i[i]])[w],c(1,2,3,the_i[i],6)]
  
  snpsOfInterest<-as.character(all_good_chrom_a[all_good_chrom_a$genes%in%genes_g$gene,"snp_id"])


  pdf(file = paste("manhattan_plot_lfmm_candidate_genes/manhattan_plots_A.aegypti_q0.01_genes",test_names_env[i], ".pdf", sep=""), width = 16)
  manhattan(all_good_chrom_a,bp ="pos",chr = "chr", p = test_names_env[i], snp = "snp_id", logp = TRUE, ylim = c(0,15), genomewideline = 4.8, suggestiveline = FALSE ,xlab = "chromosome",highlight = snpsOfInterest,
            ylab = "-log10(P)", cex= 1.1, cex.axis=1.5, cex.lab=1.3,col = c("blue4","orange3"),main = test_names_env[i])     
  dev.off()
}


chromosome<-unique(all_good_chrom_lfmm$chr)
for (i in 1:length (the_i)){
  
  for (j in 1:length(chromosome)){
  
  genes_g<-genes[genes$var==test_names_env[i] & genes$CHR==j,]
  genes_zoom<-data.frame("gene"=genes_g[genes_g$CHR==j,5])
  genes_g<-na.omit(genes_g)
  
  all_good_chrom_b<-all_good_chrom_lfmm[all_good_chrom_lfmm$chr==chromosome[j],c(1,2,3,the_i[i],29)]
  
  
  #L <- as.numeric(length(all_good_chrom_a[,the_i[i]]))
  #w = which(sort(all_good_chrom_a[,the_i[i]]) < q * (1:L) / L)
  #candidates = all_good_chrom_lfmm[order(all_good_chrom_lfmm[,the_i[i]])[w],c(1,2,3,the_i[i],6)]
  
  snpsOfInterest<-as.character(all_good_chrom_b[all_good_chrom_b$genes%in%genes_g$gene,"snp_id"])
  
  
  pdf(file = paste("manhattan_plot_lfmm_candidate_genes/manhattan_plots_lfmm_FDR0.01_A.aegypti_all_whole_genome_chrom",chromosome[j],"_topgenes_",test_names_env[i], ".pdf", sep=""), width = 16)
  manhattan(subset(all_good_chrom_b,chr == chromosome[j]),bp ="pos",chr = "chr", p = test_names_env[i], snp = "snp_id", logp = TRUE, ylim = c(0,15), genomewideline = 4.8, suggestiveline = FALSE ,xlab = "chromosome",highlight = snpsOfInterest,
            ylab = "-log10(P)", cex= 1.1, cex.axis=1.5, cex.lab=1.3,col = c("blue4","orange3"),main = test_names_env[i])     
  dev.off()
}
