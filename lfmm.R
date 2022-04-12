
#!/usr/bin/env Rscript
#rm(list = ls())
setwd("~/Documents/malaria_vector_reserach/input_data/a.aegypti/lfmm")
## genomic inflation factor: http://rstudio-pubs-static.s3.amazonaws.com/9743_8a5f7ba3aa724d4b8270c621fdf6d06e.html

var_list<-c("latitude_e","longitude_e","MAT_e","MWMT_e","MCMT_e","TD_e","MAP_e","MSP_e","AHM_e","SHM_e","DD_0_e","DD5_e","DD_18_e","DD18_e","NFFD_e","bFFP_e","eFFP_e","FFP_e","PAS_e","EMT_e","EXT_e","Eref_e","CMD_e","MAR_e","RH_e")


################################
### loop through each covariate
for (j in 22:25){ 
    z.table = NULL  
    for (i in 1:3){ ### loop through each run
    file.name = paste("cov",j,"_run_1-3/lfmm_RUN",i,"_k3_all_loci_s",j,".3.zscore", sep="")
    z.table = cbind(z.table, read.table(file.name)[,1])  ## combining zscore for each run
    }
    z.score = apply(z.table, MARGIN = 1, median) #combines z-scores
    lambda = median(z.score^2)/0.456   ## [1] 1.066392   ## round2 [1] 1.064335

    adjusted.p.values = data.frame(pchisq(z.score^2/lambda, df = 1, lower = F))    #re-adjust p-values
    colnames(adjusted.p.values)<-var_list[j]
    write.table(adjusted.p.values,file=paste("A.aegypti_K3_adjusted.p.values.cov",j,"_runs_1-3.txt", sep = ""),col.names= TRUE,row.names = FALSE, quote= FALSE)

    jpeg(file= paste("qq_plot_Pvalues_hist_lfmm_cov",j,"_runs_1-3.jpeg", sep = ""), height= 500, width = 1200)
    par(mfrow=c(1,2))
    aa<-adjusted.p.values[,1]
    qqplot(rexp(length(aa), rate = log(10)),
       -log10(aa), xlab = "Expected quantile",
       pch = 19, cex = .4, main = "qqplot of adjusted P-values")
       abline(0,1)

    hist(adjusted.p.values[,1],xlab ="adjusted p-values",main = "histogram of adjusted p-values")
    dev.off()

} ### close jloop









## cov1 (latitude)
z.table = NULL
for (i in 1:3){
  file.name = paste("cov1_run_1-3/lfmm_RUN",i,"_k3_all_loci_s1.3.zscore", sep="")
  z.table = cbind(z.table, read.table(file.name)[,1])  ## combining zscore for each run
}
z.score = apply(z.table, MARGIN = 1, median) #combines z-scores
lambda = median(z.score^2)/0.456   ## [1] 1.066392   ## round2 [1] 1.064335

adjusted.p.values.cov1 = data.frame(pchisq(z.score^2/lambda, df = 1, lower = F))    #re-adjust p-values

colnames(adjusted.p.values.cov1)<-"latitude_e"
write.table(adjusted.p.values.cov1,file="A.aegypti_K3_adjusted.p.values.cov1_runs_1-3.txt",col.names= TRUE,row.names = FALSE, quote= FALSE)

hist(adjusted.p.values.cov1$latitude_e,col="red",xlab ="adjusted p-values",main = "histogram of adjusted p-values-latitude")

### look at the qq-plot and distribution of pvalues
jpeg(file= "qq_Pvalues_hist_lfmm_latitude_runs_1-3.jpeg", height= 500, width = 1200)
par(mfrow=c(1,2))
aa<-adjusted.p.values.cov1$latitude_e
qqplot(rexp(length(aa), rate = log(10)),
       -log10(aa), xlab = "Expected quantile",
       pch = 19, cex = .4)
abline(0,1)

hist(adjusted.p.values.cov1$latitude_e,col="red",xlab ="adjusted p-values",main = "histogram of adjusted p-values-latitude")

dev.off()





#####
#!/usr/bin/env Rscript
library("qvalue")
setwd("/share/lanzarolab/users/shaghayegh/data/A.aegypti/good_ind_Q30.DP5.caling_1_maf0.03/downstream_analysis/lfmm")
all_good_chrom_lfmm<-read.table(file = "var_out_A.aegypti_Q30.DP5.caling_1_maf0.03_lfmm_IS_STD_12_pops.txt", header = FALSE)
colnames(all_good_chrom_lfmm)<-c("snp_id","chr","pos","latitude_e","longitude_e","MAT_e","MWMT_e","MCMT_e","TD_e","MAP_e","MSP_e","AHM_e","SHM_e","DD_0_e","DD5_e","DD_18_e","DD18_e","NFFD_e","bFFP_e","eFFP_e","FFP_e","PAS_e","EMT_e","EXT_e","Eref_e","CMD_e","MAR_e","RH_e")


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
  candidates = all_good_chrom_lfmm[order(all_good_chrom_lfmm[,the_i[i]])[w],c(1,2,3,the_i[i])]
  write.table(candidates, file = paste("lfmm_A.aegypti_12_pops_list_candidate_loci_FDR_0.01_",test_names_env[i],".txt", sep = ""))
}


### gives the same results as above
## q-values
#aa$qval <- qvalue(aa$latitude_e)$qvalues
#alpha <- 0.01
#outliers <- aa[which(aa$qval < alpha),]


## Benjamin-Hochberg (similar to LFMM)
#pc_data$padj <- p.adjust(pc_data$Pvalue,method="BH")
#alpha <- 0.01
#outliers_BH <- pc_data[which(padj < alpha),]






