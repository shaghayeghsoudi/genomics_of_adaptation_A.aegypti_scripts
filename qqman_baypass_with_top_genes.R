#!/usr/bin/env Rscript

rm(list = ls())
library(qqman)

setwd("/share/lanzarolab/users/shaghayegh/data/A.aegypti/good_ind_Q30.DP5.caling_1_maf0.03/downstream_analysis/baypass/iS_STD_covariate_model_outlier_analysis")
#setwd("~/Documents/malaria_vector_reserach/input_data/a.aegypti/baypass")


all_good_chrom<-read.table(file = "var_out_A.aegypti_Q30.DP5.caling_1_maf0.03_BayPass_IS_STD_12_pops_with_genes_attached.txt",header = FALSE)
#all_good_chrom<-read.table(file = "var_out_test.txt",header = FALSE)


colnames(all_good_chrom)<-c("snp_id","CHR","pos","latitude_e","longitude_e","MAT_e","MWMT_e","MCMT_e","TD_e","MAP_e","MSP_e","AHM_e","SHM_e","DD_0_e","DD5_e","DD_18_e","DD18_e","NFFD_e","bFFP_e","eFFP_e","FFP_e","PAS_e","EMT_e","EXT_e","Eref_e","CMD_e","MAR_e","RH_e","genes")
#colnames(all_good_chrom)<-c("snp_id","CHR","pos","latitude_e","longitude_e","genes")


test_type_env<-array ("0",ncol (all_good_chrom))
env_type<-grep ("_e$",colnames(all_good_chrom))
test_type_env[min(env_type):(min(env_type)+24)]<-"envir"
all_good_chrom$id<- paste (all_good_chrom[,2],all_good_chrom[,3], sep = "__")

#### Top candidate windows and BF threshold ####
out_res <- NULL
test_names_env <- colnames (all_good_chrom)[test_type_env == "envir"]
the_i <- which (test_type_env == "envir") 
bf <- 20


## extract outliers (here SNPs with BFs more than significance threshold based on Baypass paper )

genes<-read.table(file = "baypass_A.aegypi_12pops_BF20_list_candidate_loci_with_attached_genes_U_ALL_Vars.table", header = TRUE)
colnames(genes)<-c("snp_id","CHR","pos","BF_is","genes","var")
chromosome<-unique(all_good_chrom$CHR)

for(i in 1:length(test_names_env)){
  
  genes_g<-data.frame("gene"=genes[genes$var==test_names_env[i],5])
  genes_g<-na.omit(genes_g)
  #genes_zoom<-data.frame("gene"=genes_g[genes_g$chrom==j,2])
  all_good_chrom_a<-all_good_chrom[,c(1,2,3,the_i[i],29)]
  snpsOfInterest<-as.character(all_good_chrom_a[all_good_chrom_a$genes%in%genes_g$gene,"snp_id"])
  
  
  ###### plot convergent genes on the manhattan plot (across all chromosomes) #####
  pdf(file = paste("manhattan_plot_baypass_outlier_genes/manhattan_plots_baypass_STD_IS2_A.aegypti_all_whole_genome_topgenes_",test_names_env[i], ".pdf", sep=""), width = 18)
  manhattan(all_good_chrom_a,bp ="pos", chr = "CHR", p = test_names_env[i], snp = "snp_id", logp = FALSE, ylim = c(0,60), genomewideline = 20, suggestiveline = FALSE ,xlab = "chromosome",highlight = snpsOfInterest,
            ylab = "BFis (in dB)", cex = 1.1, cex.axis=1.5, cex.lab=1.4,col = c("blue4","orange3"),main = test_names_env[i])     
  dev.off()
  
}


###


for(i in 1:length(test_names_env)){
  
  for (j in 1:length(chromosome)){
  
  
  genes_g<-genes[genes$var==test_names_env[i] & genes$CHR==j,]
  genes_zoom<-data.frame("gene"=genes_g[genes_g$CHR==j,5])
  genes_zoom<-na.omit(genes_zoom)
  all_good_chrom_b<-all_good_chrom[all_good_chrom$CHR==chromosome[j],c(1,2,3,the_i[i],6)]
  snpsOfInterest<-as.character(all_good_chrom_b[all_good_chrom_b$genes%in%genes_zoom$gene,"snp_id"])
  
  
  ###### plot convergent genes on the manhattan plot (across all chromosomes) #####
  pdf(file = paste("manhattan_plot_baypass_outlier_genes/manhattan_plots_baypass_STD_IS2_A.aegypti_all_whole_genome_chrom",chromosome[j],"_topgenes_",test_names_env[i], ".pdf", sep=""), width = 16)
  manhattan(subset(all_good_chrom_b,CHR == chromosome[j]), bp ="pos", p = test_names_env[i], snp = "snp_id", logp = FALSE, ylim = c(0,60), genomewideline = 20, suggestiveline = FALSE ,xlab = "chromosome",highlight = snpsOfInterest,
            ylab = "BFis (in dB)", cex = 1.1, cex.axis=1.5, cex.lab=1.4,col = c("blue4","orange3"),main = test_names_env[i])     
  dev.off()
  
  }
}
