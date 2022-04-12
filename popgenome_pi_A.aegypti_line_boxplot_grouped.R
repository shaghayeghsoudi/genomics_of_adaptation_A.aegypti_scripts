#https://github.com/parkingvarsson/PhotoperiodLocalAdaptation/blob/master/6-Positive%20selection/chr10_genome.sig.angsd_tP_tajD.group_plot.R

setwd("~/Documents/malaria_vector_reserach/input_data/a.aegypti/baypass")

mvmt<-read.table(file = "baypass_A.aegypi_12pops_BF20_list_candidate_loci_with_attached_genes_U_MWMT_e.txt", header = TRUE)

mvmt_table<-data.frame (table (as.character (mvmt[,"genes"])))
mvmt_table10<-mvmt_table[mvmt_table$Freq>=10,]
#aa<-mvmt_table10$Var1

### load genes bed file
g_coor<-read.table(file = "genes_good.bed.txt", header = FALSE)
colnames(g_coor)<-c("chrom","start","stop","gene")

#g_coor$start<-g_coor$start/10^6
#g_coor$stop<-g_coor$stop/10^6
#g_coor_gene<-unique(g_coor$gene)

genes_top10<-g_coor[g_coor$gene%in%mvmt_table10$Var1,]
g_coor_gene<-genes_top10$gene


### load diversity_neutrality data
ndch1<-read.table(file = "diversity_neutrality_grouped_3_pops/A.aegypti_nd_chrom1_12pops_3Groups_50k_5K.txt", header = TRUE)
ndch1$chrom<-rep("1",nrow(ndch1))

ndch2<-read.table(file = "diversity_neutrality_grouped_3_pops/A.aegypti_nd_chrom2_12pops_3Groups_50k_5K.txt", header = TRUE)
ndch2$chrom<-rep("2",nrow(ndch2))

ndch3<-read.table(file = "diversity_neutrality_grouped_3_pops/A.aegypti_nd_chrom3_12pops_3Groups_50k_5K.txt", header = TRUE)
ndch3$chrom<-rep("3",nrow(ndch3))

data_pre<-do.call("rbind", list(ndch1, ndch2, ndch3))
data_pre<-na.omit(data_pre)


### find thresholds ####
### quantile (pi)
the_quantile <- 0.05
quantile_pi_G1 <- quantile (data_pre[,"A.ae1_pi"], the_quantile,na.rm = T)
quantile_pi_G2 <- quantile (data_pre[,"A.ae2_pi"], the_quantile,na.rm = T)
quantile_pi_G3 <- quantile (data_pre[,"A.ae3_pi"], the_quantile,na.rm = T)




pdf(file = "A.aegypti_grouped3_pops_Pi.pdf", width = 8, height = 4)
for (i in 1:length(g_coor_gene)){
  
  ## find the gene in stats file
  genes_top10_focal<-genes_top10[genes_top10$gene==g_coor_gene[i],]
  aa<-data_pre[data_pre$chrom==genes_top10_focal$chrom,]
  bb<-aa[aa$start>genes_top10_focal$start & aa$stop<genes_top10_focal$stop,]
  
  
  ### select the region for plot
  region<-aa[aa$start>(genes_top10_focal$start)-1000000 & aa$stop<(genes_top10_focal$stop)+1000000,]
  #region$mid_good<-region$mid/1000000
 
  
  plot(region$mid,region$A.ae1_pi , type = "l",lwd = 1.8, col = "darkcyan", xlab = "Position (Mb)",las = 1,ylab="",  ylim=c(0,0.002), main=paste("A.aegypti_",g_coor_gene[i]))
  mtext(side=2, text="pi", line=4.5)
  
  lines(region$mid, region$A.ae2_pi, type="l",lwd = 1.5, col = "dodgerblue4")
  lines(region$mid, region$A.ae3_pi, type="l",lwd = 1.5, col = "sienna3")
  
  
  
  abline(v =c(genes_top10_focal$start/1000000,genes_top10_focal$stop/1000000) ,col="red", lwd=1.4, lty=2)
  legend("topleft",legend = c("A.ae_G1","A.ae_G2","A.ae_G3"),col = c("darkcyan","dodgerblue4","sienna3"), lwd=0.8, cex = 0.8, bg = "white",bty="n",lty = 1,horiz=FALSE,y.intersp=0.9) 
  #abline(v = c(g_coor$start, g_coor$stop), col="orangered4", lwd=2, lty=2)
  abline(h =c(quantile_pi_G1,quantile_pi_G2,quantile_pi_G3) ,col=c("darkcyan","dodgerblue4","sienna3"), lwd=0.5, lty=4)
  
}
dev.off()


##################
### box plots (pi) ####


pdf(file = "A.aegypti_grouped3_pops_Pi_Boxplots_true.pdf", width = 7, height = 5)
for (i in 1:length(g_coor_gene)){
  
  ## find the gene in stats file
  genes_top10_focal<-genes_top10[genes_top10$gene==g_coor_gene[i],]
  aa<-data_pre[data_pre$chrom==genes_top10_focal$chrom,]
  bb<-aa[aa$start>genes_top10_focal$start & aa$stop<genes_top10_focal$stop,]
  
  
  ### select the region for plot
  region<-aa[aa$start>(genes_top10_focal$start) & aa$stop<(genes_top10_focal$stop),]
  nonregion<-aa[!(aa$start>(genes_top10_focal$start) & aa$stop<(genes_top10_focal$stop)),]
  #region$mid_good<-region$mid/1000000
  
  ## G1
  region_G1<-region[,c("start","stop","mid","A.ae1_pi","chrom")]
  region_G1$group<-rep("A.ae1",nrow(region_G1))
  region_G1$status<-rep("wholegenome",nrow(region_G1))
  
  nonregion_G1<-nonregion[,c("start","stop","mid","A.ae1_pi","chrom")]
  nonregion_G1$group<-rep("A.ae1",nrow(nonregion_G1))
  nonregion_G1$status<-rep("focal",nrow(nonregion_G1))
  G1<-do.call("rbind", list(region_G1,nonregion_G1))
  colnames(G1)<-c("start","stop","mid", "A.ae_pi","chrom","group","status")
  
  ## G2
  region_G2<-region[,c("start","stop","mid","A.ae2_pi","chrom")]
  region_G2$group<-rep("A.ae2",nrow(region_G2))
  region_G2$status<-rep("wholegenome",nrow(region_G2))
  
  nonregion_G2<-nonregion[,c("start","stop","mid","A.ae2_pi","chrom")]
  nonregion_G2$group<-rep("A.ae2",nrow(nonregion_G2))
  nonregion_G2$status<-rep("focal",nrow(nonregion_G2))
  G2<-do.call("rbind", list(region_G2,nonregion_G2))
  colnames(G2)<-c("start","stop","mid", "A.ae_pi","chrom","group","status")
  
 
   ## G3
  region_G3<-region[,c("start","stop","mid","A.ae3_pi","chrom")]
  region_G3$group<-rep("A.ae3",nrow(region_G3))
  region_G3$status<-rep("wholegenome",nrow(region_G3))
  
  nonregion_G3<-nonregion[,c("start","stop","mid","A.ae3_pi","chrom")]
  nonregion_G3$group<-rep("A.ae3",nrow(nonregion_G3))
  nonregion_G3$status<-rep("focal",nrow(nonregion_G3))
  G3<-do.call("rbind", list(region_G3,nonregion_G3))
  colnames(G3)<-c("start","stop","mid", "A.ae_pi","chrom","group","status")
 
  allGroup <- do.call("rbind", list(G1,G2,G3)) 
  
  
  
  
  ###1.3.4 Making the plot 
  #diversity
  #lvl1=c("chr10-700kb","genome-wide")
  lvl2=c("A.aeG1","A.aeG2","A.aeG3")
  #factor1=as.factor(c(rep("chr10-700kb",117),rep("genome-wide",38003),rep("chr10-700kb",117),rep("genome-wide",38003),rep("chr10-700kb",117),rep("genome-wide",38001)))
  #factor2=as.factor(c(rep("North",38120),rep("Mid",38120),rep("South",38118)))
 # plotgrp=factor(paste(factor2,factor1),levels=c(sapply(lvl2,paste,lvl1)))
  
  
  
  cols=c("darkcyan","dodgerblue4","sienna3")
  at_1=c(1:2,4:5,7:8)
  at=c(1.5,4.5,7.5)
  
  
  par(las=1)
  par(mar=c(3,8,2,1))
          
  
  lvl2=c("A.aeG1","A.aeG2","A.aeG3")
  boxplot(A.ae_pi~group*status, data=allGroup,notch=FALSE,at=at_1,xaxt="n",cex.axis=1.5,las = 3,frame.plot = FALSE,xlab="",ylab="",cex.lab=2.5,labs=2,col=c(cols[3],"grey20",cols[2],"grey50",cols[1],"grey80"),outline=FALSE,ylim=c(0,0.002), main =paste("A.aegypti_",g_coor_gene[i]) )
  axis(1,at=at,labels=lvl2,tick=T,cex.axis=1.5)
  mtext(expression(pi),las=1,side=2,line=4,cex=3)
  #legend("topleft", fill = rainbow(3, s = 0.5), legend = c(1,2,3), horiz = T)
  segments(1,0.032,1,0.033)
  segments(1,0.033,1.95,0.033)
  segments(1.95,0.032,1.95,0.033)
  text(1.5,0.034,"***",cex=2)
  segments(4,0.032,4,0.033)
  segments(4,0.033,5,0.033)
  segments(5,0.032,5,0.033)
  text(4.5,0.034,"***",cex=2)
  segments(7,0.032,7,0.033)
  segments(7,0.033,8,0.033)
  segments(8,0.032,8,0.033)
  text(7.5,0.034,"***",cex=2)
  
  
  
}
dev.off()

### box plots


