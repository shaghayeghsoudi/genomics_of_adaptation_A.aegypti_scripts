#!/usr/bin/env Rscript
#rm(list = ls())

library(qqman)
#######################################
##### generate Fst manhattan plots ######
#######################################

setwd("~/Documents/malaria_vector_reserach/input_data/introgression_mali_benin/Fst/windowed_nonoverlaping_50K")

## no window Fst estimations
chrom_2L<-read.table(file = "2L_Benin_Col_Gam_no_window_out.weir.fst", header = TRUE)
chrom_2R<-read.table(file = "2R_Benin_Col_Gam_no_window_out.weir.fst", header = TRUE)
chrom_3L<-read.table(file = "3L_Benin_Col_Gam_no_window_out.weir.fst", header = TRUE)
chrom_3R<-read.table(file = "3R_Benin_Col_Gam_no_window_out.weir.fst", header = TRUE)
chrom_x<-read.table(file = "X_Benin_Col_Gam_no_window_out.weir.fst", header = TRUE)

## combine all chroms
all_chrom <- do.call("rbind", list(chrom_2L, chrom_2R, chrom_3L,chrom_3R,chrom_x)) 

all_chrom$snp_id<-paste(all_chrom$CHROM,all_chrom$POS,sep ="__")
all_chrom<-na.omit(all_chrom)

all_chrom<-all_chrom[,c(4,1,2,3)]
colnames(all_chrom)<-c("SNP","CHR","BP","P")

levels(all_chrom$CHR) <- sub("2L", "1", levels(all_chrom$CHR))
levels(all_chrom$CHR) <- sub("2R", "2", levels(all_chrom$CHR))
levels(all_chrom$CHR) <- sub("3L", "3", levels(all_chrom$CHR))
levels(all_chrom$CHR) <- sub("3R", "4", levels(all_chrom$CHR))
levels(all_chrom$CHR) <- sub("X", "5", levels(all_chrom$CHR))
all_chrom$CHR<-as.numeric(as.character(all_chrom$CHR))

pdf(file = paste("manhattan_plots_Fst_Benin_An.gambiae_An.coluzzii_no_window.pdf", sep=""), width = 8)
manhattan(all_chrom,  ylim = c(0,1.1), logp = FALSE,genomewideline = F, suggestiveline = FALSE ,xlab = "chromosome",
          ylab = "Fst", cex.axis=1.5, cex.lab=1.4,col = c("blue4","orange3"),main = "Benin_A.coluzzii_A.gambiae",chrlabs =c("2L","2R","3L","3R","X")) 

dev.off()





################################################################################################
################################################################################################

## manhattan plot with ggplot2
## https://www.r-graph-gallery.com/101_Manhattan_plot.html
library("dplyr")
library("ggplot2")

colnames(chrom)<-c("SNP","CHR","BP","P")

#First of all, we need to compute the cumulative position of SNP.

don <- chrom %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(chrom, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot)


#Then we need to prepare the X axis. Indeed we do not want to display the cumulative position of SNP in bp, but just show the chromosome name instead.
axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


#Ready to make the plot using ggplot2:
ggplot(don, aes(x=BPcum, y=P)) +
  
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
















https://danielroelfs.com/blog/how-i-create-manhattan-plots-using-ggplot/
  manhplot <- ggplot(gwas.dat, aes(x = BPcum, y = -log10(P), 
                                   color = as.factor(CHR), size = -log10(P))) +
  geom_point(alpha = 0.75) +
  geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "dashed") + 
  scale_x_continuous(label = axis.set$CHR, breaks = axis.set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
  scale_color_manual(values = rep(c("#276FBF", "#183059"), nCHR)) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL, 
       y = "-log10(p)") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5)
  )

#########################################
#########################################
## line plot Fst ##
rm(list = ls())
setwd("~/Documents/malaria_vector_reserach/input_data/introgression_mali_benin/data_March_2020/1-3MB_2L_phased_50K_10Kstepping")
aa<-read.table(file = "2L_3MB_Benin_Col_Gam_out.windowed.weir.fst", header = TRUE)
#plot(aa$MEAN_FST~aa$BIN_END,type="l",col="red")

bb<-read.table(file ="2L_3MB_Mali_Col_Gam_out.windowed.weir.fst", header = TRUE)
#plot(bb$MEAN_FST~bb$BIN_END,type="l",col="red")


###### plot MEAN Fst ########
### find 0.99 quantile
the_quantile_td_p <- 0.99
quantile_aa <- quantile (aa[,"MEAN_FST"], the_quantile_td_p,na.rm = T)
quantile_bb <- quantile (bb[,"MEAN_FST"], the_quantile_td_p,na.rm = T)

start<-(2358158-250)
end<-(2431617+250)

### plot them together
pdf(file = "2L_3MB_An.gambiae_An.coluzzii_Benin_Mali_Mean_Weighted_Fst_with_kdr.pdf", width = 9, height= 8)
par(mfrow=c(2,1))
plot(aa$MEAN_FST~aa$BIN_END, type = "l",lwd = 2, col = "forestgreen", xlab = "Position ",cex.lab= 1.3,las = 1,ylab="Fst",  ylim=c(0,1),cex.axis=1.3,main = "2L_3MB_An.gambiae vs. An.coluzzii(Mean)")
mtext(side=2, text="Fst", line=5)         

lines(bb$MEAN_FST~bb$BIN_END, type="l",lwd = 1.8, col = "darkred")
legend("topleft",legend = c("Benin","Mali"),col = c("forestgreen","darkred"), cex = 1.1, bg = "white",bty="n",lty = 1:1,horiz=FALSE,x.intersp=0.8)
abline(h =c(quantile_aa,quantile_bb) ,col=c("forestgreen","darkred"), lwd=1.4, lty=2:3)
abline(v =c(start,end) ,col="red", lwd=1.4, lty=2)

#dev.off()

###### plot MEAN Fst ########

the_quantile_td_p <- 0.99
quantile_aa <- quantile (aa[,"WEIGHTED_FST"], the_quantile_td_p,na.rm = T)
quantile_bb <- quantile (bb[,"WEIGHTED_FST"], the_quantile_td_p,na.rm = T)

### plot them together
#pdf(file = "2L_3MB_An.gambiae_An.coluzzii_Benin_Mali_WEIGHTED_FST.pdf", width = 9, height= 4)
plot(aa$WEIGHTED_FST~aa$BIN_END, type = "l",lwd = 2, col = "forestgreen", xlab = "Position",cex.lab= 1.3,las = 1,ylab="Fst",  ylim=c(0,1),cex.axis=1.3,main = "2L_3MB_An.gambiae vs. An.coluzzii (weighted)")
mtext(side=2, text="Fst", line=5)         

lines(bb$WEIGHTED_FST~bb$BIN_END, type="l",lwd = 1.8, col = "darkred")
legend("topleft",legend = c("Benin","Mali"),col = c("forestgreen","darkred"), cex = 1.1, bg = "white",bty="n",lty = 1:1,horiz=FALSE,x.intersp=0.8)
abline(h =c(quantile_aa,quantile_bb) ,col=c("forestgreen","darkred"), lwd=1.4, lty=2:3)
abline(v =c(start,end) ,col="red", lwd=1.4, lty=2)

dev.off()


##############################################
################# plot pi ####################
##############################################
#rm(list = ls())

setwd("~/Documents/malaria_vector_reserach/input_data/introgression_mali_benin/Fst/windowed_nonoverlaping_50K")
aa<-read.table(file = "Benin_coluzzii_windowed_50K_pi", header = TRUE)
#plot(aa$PI~aa$BIN_END,type="l",col="red",ylim=c(0,0.005))

bb<-read.table(file = "Mali_coluzzii_windowed_50K_pi", header = TRUE)
#plot(bb$PI~bb$BIN_END,type="l",col="red",ylim=c(0,0.005))

cc<-read.table(file = "Benin_gambiae_windowed_50K_pi", header = TRUE)

dd<-read.table(file = "Mali_gambiae.windowed_50K.pi", header = TRUE)


#
the_quantile_td_p <- 0.99
quantile_aa <- quantile (aa[,"PI"], the_quantile_td_p,na.rm = T)
quantile_bb <- quantile (bb[,"PI"], the_quantile_td_p,na.rm = T)
quantile_cc <- quantile (cc[,"PI"], the_quantile_td_p,na.rm = T)
quantile_dd <- quantile (dd[,"PI"], the_quantile_td_p,na.rm = T)


### plot them together
pdf(file = "2L_An.gambiae_Benin vs. An.gambiae_Mali_Pi.pdf", width = 9, height= 4)
plot(aa$PI~aa$BIN_END, type = "l",lwd = 1.5, col = "dodgerblue4", xlab = "Position (Mb)",cex.lab= 1.3,las = 1,ylab="nucleotide diversity",  ylim=c(0,0.005),cex.axis=1,main = "2L An.coluzzii (Benin)")
mtext(side=2, text="Fst", line=5)         

lines(bb$PI~bb$BIN_END, type="l",lwd = 1.2, col = "mediumvioletred")
legend("topleft",legend = c("An.gambiae (Benin)","An.gambiae (Mali)"),col = c("dodgerblue4","mediumvioletred"), cex = 0.9, bg = "white",bty="n",lty = 1:1,horiz=FALSE,x.intersp=0.8)
abline(h =c(quantile_aa,quantile_bb) ,col=c("dodgerblue4","mediumvioletred"), lwd=1.4, lty=2:3)
abline(v =c(start,end) ,col="red", lwd=1.4, lty=2)

dev.off()




### plot nucleotide dicersity seperately
pdf(file = "2L_Pi_each_species_seperately.pdf", width = 9, height= 12)

par(mfrow=c(4,1))
plot(aa$PI~aa$BIN_END, type = "l",lwd = 1.5, col = "dodgerblue4", xlab = "Position",cex.lab= 1.3,las = 1,ylab="nucleotide diversity",  ylim=c(0,0.005),cex.axis=1,main = "2L An.coluzzii (Benin)")
mtext(side=2, text="Fst", line=5)         
abline(h =c(quantile_aa) ,col=c("dodgerblue4"), lwd=1.4, lty=2:3)
abline(v =c(start,end) ,col="red", lwd=1.4, lty=2)


plot(bb$PI~bb$BIN_END, type = "l",lwd = 1.5, col = "chartreuse4", xlab = "Position",cex.lab= 1.3,las = 1,ylab="nucleotide diversity",  ylim=c(0,0.005),cex.axis=1,main = "2L An.coluzzii (Mali)")
mtext(side=2, text="Fst", line=5)         
abline(h =c(quantile_bb) ,col=c("chartreuse4"), lwd=1.4, lty=2:3)
abline(v =c(start,end) ,col="red", lwd=1.4, lty=2)


plot(cc$PI~cc$BIN_END, type = "l",lwd = 1.5, col = "darkorange3", xlab = "Position",cex.lab= 1.3,las = 1,ylab="nucleotide diversity",  ylim=c(0,0.005),cex.axis=1,main = "2L An.gambiae (Benin)")
mtext(side=2, text="Fst", line=5)         
abline(h =c(quantile_cc) ,col=c("darkorange3"), lwd=1.4, lty=2:3)
abline(v =c(start,end) ,col="red", lwd=1.4, lty=2)


plot(dd$PI~dd$BIN_END, type = "l",lwd = 1.5, col = "deeppink4", xlab = "Position",cex.lab= 1.3,las = 1,ylab="nucleotide diversity",  ylim=c(0,0.005),cex.axis=1,main = "2L An.gambiae (Mali)")
mtext(side=2, text="Fst", line=5)         
abline(h =c(quantile_dd) ,col=c("deeppink4"), lwd=1.4, lty=2:3)
abline(v =c(start,end) ,col="red", lwd=1.4, lty=2)

dev.off()



#####  plots nucleotide diversity (pi) in a smoothing form plot ####
#rm(list = ls())
setwd("~/Documents/malaria_vector_reserach/input_data/introgression_mali_benin/pi/3R")

data_a<-read.table(file = "Mali_gambiae_3R_50K_10K.windowed.pi", header = TRUE)
data_b<-read.table(file = "Mali_gambiae_3R_50K_100.windowed.pi", header = TRUE)

pdf(file = "3R_Pi_smoothing_An.gambiae_Mali.pdf", width = 8, height= 8)
plot(PI~BIN_END,data=data_a,ylab="Nucleotide diversity", xlab= "Position", cex.axis=1.1,main = "3R An.gambiae (Mali)",cex.lab= 1.3,ylim=c(0,0.005) )
with(data_b, lines(lowess(BIN_END,PI), col="tomato", lwd=4))
dev.off()


data_a<-read.table(file = "Mali_coluzzii_3R_50K_10K.windowed.pi", header = TRUE)
data_b<-read.table(file = "Mali_coluzzii_3R_50K_100.windowed.pi", header = TRUE)

pdf(file = "3R_Pi_smoothing_An.coluzzii_Mali.pdf", width = 8, height= 8)
plot(PI~BIN_END,data=data_a,ylab="Nucleotide diversity", xlab= "Position", cex.axis=1.1,ylim=c(0,0.005),main = "3R An.coluzzii (Mali)",cex.lab= 1.3 )
with(data_b, lines(lowess(BIN_END,PI), col="tomato", lwd=4))
dev.off()


data_a<-read.table(file = "Benin_coluzzii_3R_50K_10K.windowed.pi", header = TRUE)
data_b<-read.table(file = "Benin_coluzzii_3R_50K_100.windowed.pi", header = TRUE)

pdf(file = "3R_Pi_smoothing_An.coluzzii_Benin.pdf", width = 8, height= 8)
plot(PI~BIN_END,data=data_a,ylab="Nucleotide diversity", xlab= "Position", cex.axis=1.1,ylim=c(0,0.005),main = "3R An.coluzzii (Benin)",cex.lab= 1.3 )
with(data_b, lines(lowess(BIN_END,PI), col="tomato", lwd=4))
dev.off()




data_a<-read.table(file = "Benin_gambiae_3R_50K_10K.windowed.pi", header = TRUE)
data_b<-read.table(file = "Benin_gambiae_3R_50K_100.windowed.pi", header = TRUE)

pdf(file = "3R_Pi_smoothing_An.gambiae_Benin.pdf", width = 8, height= 8)
plot(PI~BIN_END,data=data_a,ylab="Nucleotide diversity", xlab= "Position", cex.axis=1.1,ylim=c(0,0.005),main = "3R An.gambiae (Benin)",cex.lab= 1.3 )
with(data_b, lines(lowess(BIN_END,PI), col="tomato", lwd=4))
dev.off()


### plots smoothing plots in comparison ###
### (A) compare Benin coluzzii and Mali coluzzii

pdf(file = "3R_Pi_smoothing_all_type_comparisons.pdf", width = 12, height= 12)
par(mfrow=c(2,2))
data_a<-read.table(file = "Benin_coluzzii_3R_50K_10K.windowed.pi", header = TRUE)
data_b<-read.table(file = "Benin_coluzzii_3R_50K_100.windowed.pi", header = TRUE)
data_c<-read.table(file = "Mali_coluzzii_3R_50K_100.windowed.pi", header = TRUE)

#pdf(file = "2L_Pi_smoothing_An.coluzzii_Benin_Benin_Mali.pdf", width = 8, height= 8)
plot(PI~BIN_END,data=data_a,ylab="Nucleotide diversity", xlab= "Position", cex.axis=1.1,ylim=c(0,0.005),main = "3R An.coluzzii(B)-An.coluzzii(M)",cex.lab= 1.3 )
with(data_b, lines(lowess(BIN_END,PI), col="tomato", lwd=4))
with(data_c, lines(lowess(BIN_END,PI), col="chartreuse4", lwd=4))
legend("topleft",legend = c("A.coluzzii (Benin)","A.coluzzii (Mali)"),col = c("tomato","chartreuse4"), lwd=1, cex = 0.9, bg = "white",bty="n",lty = 1,horiz=FALSE,y.intersp=0.9) 
#dev.off()

####
### (B) compare Benin gambiae and Mali gambiae
data_a<-read.table(file = "Benin_gambiae_3R_50K_10K.windowed.pi", header = TRUE)
data_b<-read.table(file = "Benin_gambiae_3R_50K_100.windowed.pi", header = TRUE)
data_c<-read.table(file = "Mali_gambiae_3R_50K_100.windowed.pi", header = TRUE)

#pdf(file = "2L_Pi_smoothing_An.gambiae_Benin_Benin_Mali.pdf", width = 8, height= 8)
plot(PI~BIN_END,data=data_a,ylab="Nucleotide diversity", xlab= "Position", cex.axis=1.1,ylim=c(0,0.005),main = "3R An.gambiae(B)-An.gambiae(M)",cex.lab= 1.3 )
with(data_b, lines(lowess(BIN_END,PI), col="tomato", lwd=4))
with(data_c, lines(lowess(BIN_END,PI), col="chartreuse4", lwd=4))
legend("topleft",legend = c("A.gambiae (Benin)","A.gambiae (Mali)"),col = c("tomato","chartreuse4"), lwd=1, cex = 0.9, bg = "white",bty="n",lty = 1,horiz=FALSE,y.intersp=0.9) 
#dev.off()

### (B) compare Benin gambiae and Benin coluzzii
#rm(list = ls())
data_a<-read.table(file = "Benin_gambiae_3R_50K_10K.windowed.pi", header = TRUE)
data_b<-read.table(file = "Benin_gambiae_3R_50K_100.windowed.pi", header = TRUE)
data_c<-read.table(file = "Benin_coluzzii_3R_50K_100.windowed.pi", header = TRUE)

#pdf(file = "2L_Pi_smoothing_An.gambiae_Benin.Benin_An.coluzzii.Benin.pdf", width = 8, height= 8)
plot(PI~BIN_END,data=data_a,ylab="Nucleotide diversity", xlab= "Position", cex.axis=1.1,ylim=c(0,0.005),main = "3R An.gambiae(B)-An.coluzzii(B)",cex.lab= 1.3 )
with(data_b, lines(lowess(BIN_END,PI), col="tomato", lwd=4))
with(data_c, lines(lowess(BIN_END,PI), col="chartreuse4", lwd=4))
legend("topleft",legend = c("A.gambiae (Benin)","A.coluzzii (Benin)"),col = c("tomato","chartreuse4"), lwd=1, cex = 0.9, bg = "white",bty="n",lty = 1,horiz=FALSE,y.intersp=0.9) 
#dev.off()


### (c) compare Benin gambiae and Benin coluzzii
#rm(list = ls())
data_a<-read.table(file = "Mali_gambiae_3R_50K_10K.windowed.pi", header = TRUE)
data_b<-read.table(file = "Mali_gambiae_3R_50K_100.windowed.pi", header = TRUE)
data_c<-read.table(file = "Mali_coluzzii_3R_50K_100.windowed.pi", header = TRUE)

#pdf(file = "2L_Pi_smoothing_An.gambiae_Mali.Mali_An.coluzzii.Mali.pdf", width = 8, height= 8)
plot(PI~BIN_END,data=data_a,ylab="Nucleotide diversity", xlab= "Position", cex.axis=1.1,ylim=c(0,0.005),main = "3R An.gambiae(M)-An.coluzzii(M)",cex.lab= 1.3 )
with(data_b, lines(lowess(BIN_END,PI), col="tomato", lwd=4))
with(data_c, lines(lowess(BIN_END,PI), col="chartreuse4", lwd=4))
legend("topleft",legend = c("A.gambiae (Mali)","A.coluzzii (Mali)"),col = c("tomato","chartreuse4"), lwd=1, cex = 0.9, bg = "white",bty="n",lty = 1,horiz=FALSE,y.intersp=0.9) 

dev.off()


#######################################################################
################### plot Tajima's D ###################################
#######################################################################

setwd("~/Documents/malaria_vector_reserach/input_data/introgression_mali_benin/TjimasD/X")
tajima<-read.table(file = "X_sparrow_Tajimas_D_50K_5K_all_pops.txt", header = TRUE)

the_quantile_td <- 0.05
qpre_td_ColB <- quantile (tajima[,"A.coluzzii_Benin_TD"], the_quantile_td,na.rm = T)
qpre_td_ColM <- quantile (tajima[,"A.coluzzii_Mali_TD"], the_quantile_td,na.rm = T)
qpre_td_GamB <- quantile (tajima[,"A.gambiae_Benin_TD"], the_quantile_td,na.rm = T)
qpre_td_GamM <- quantile (tajima[,"A.gambiae_Mali_TD"], the_quantile_td,na.rm = T)

the_quantile_td_p <- 0.95
qpre2_td_ColB <- quantile (tajima[,"A.coluzzii_Benin_TD"], the_quantile_td_p,na.rm = T)
qpre2_td_ColM <- quantile (tajima[,"A.coluzzii_Mali_TD"], the_quantile_td_p,na.rm = T)
qpre2_td_GamB <- quantile (tajima[,"A.gambiae_Benin_TD"], the_quantile_td,na.rm = T)
qpre2_td_GamM <- quantile (tajima[,"A.gambiae_Mali_TD"], the_quantile_td,na.rm = T)



###plot 
pdf(file = "X_An.gambiae_An.gambiae_Tajima_D.pdf", width = 12, height= 6)
plot(tajima$stop ,tajima$A.coluzzii_Benin_TD, type = "l",lwd = 1.2, col = "darkcyan", xlab = "Position",las = 1,ylab="TD",  cex.axis=1.4,cex.lab = 1.3,ylim=c(-3,5),main = "Tajima's D - XL")
mtext(side=2, text="Tajima's D", line=4.5)

lines(tajima$stop, tajima$A.coluzzii_Mali_TD, type="l",lwd = 1.5, col = "dodgerblue4")
lines(tajima$stop, tajima$A.gambiae_Benin_TD, type="l",lwd = 1.5, col = "sienna3")
lines(tajima$stop, tajima$A.gambiae_Mali_TD, type="l",lwd = 1.5, col = "plum4")


legend("topright",legend = c("A.coluzzii(Benin)","A.coluzzii(Mali)","A.gambiae(Benin)","A.gambiae(Mali)"),col = c("darkcyan","dodgerblue4","sienna3","plum4"), lwd=0.8, cex = 0.8, bg = "white",bty="n",lty = 1,horiz=FALSE,y.intersp=0.9) 
#abline(h =c(0,qpre_td_ColB,qpre_td_ColM,qpre_td_GamB,qpre_td_GamM,qpre2_td_ColB,qpre2_td_ColM,qpre2_td_GamB,qpre2_td_GamM) ,col=c("black","darkcyan","dodgerblue4","sienna3","plum4","darkcyan","dodgerblue4","sienna3","plum4"), lwd=0.5, lty=4)
abline(h =c(0,qpre_td_ColB,qpre_td_ColM,qpre_td_GamB,qpre_td_GamM) ,col=c("black","darkcyan","dodgerblue4","sienna3","plum4"), lwd=0.7, lty=4)

dev.off()


#### plot Tajima's D seperately one plot per species (3MB) ####
setwd("~/Documents/malaria_vector_reserach/input_data/introgression_mali_benin/TjimasD/X")
tajima<-read.table(file = "X_sparrow_Tajimas_D_50K_5K_all_pops.txt", header = TRUE)


the_quantile_td <- 0.05
qpre_td_ColB <- quantile (tajima[,"A.coluzzii_Benin_TD"], the_quantile_td,na.rm = T)
qpre_td_ColM <- quantile (tajima[,"A.coluzzii_Mali_TD"], the_quantile_td,na.rm = T)
qpre_td_GamB <- quantile (tajima[,"A.gambiae_Benin_TD"], the_quantile_td,na.rm = T)
qpre_td_GamM <- quantile (tajima[,"A.gambiae_Mali_TD"], the_quantile_td,na.rm = T)


#start<-(2358158-250)
#end<-(2431617+250)


#the_quantile_td_p <- 0.95
#qpre2_td_ColB <- quantile (tajima[,"A.coluzzii_Benin_TD"], the_quantile_td_p,na.rm = T)
#qpre2_td_ColM <- quantile (tajima[,"A.coluzzii_Mali_TD"], the_quantile_td_p,na.rm = T)
#qpre2_td_GamB <- quantile (tajima[,"A.gambiae_Benin_TD"], the_quantile_td,na.rm = T)
#qpre2_td_GamM <- quantile (tajima[,"A.gambiae_Mali_TD"], the_quantile_td,na.rm = T)

pdf(file = "X_3MB_An.gambiae_An.gambiae_Tajima_D_sep.pdf", width = 10, height= 14)
par(mfrow=c(4,1))
plot(tajima$A.coluzzii_Benin_TD~tajima$stop,type="l",col="darkcyan",xlab = "Position", cex.axis=1.6,cex.lab = 1.3,ylab = "Tajima's D", main = "A.coluzzii(Benin)")
abline(h =c(0,qpre_td_ColB) ,col=c("black","darkcyan"), lwd=1.5, lty=1:4)
#abline(v = c(0,24000000), col="red", lwd=3, lty=2)
#abline(v = c(start,end), col="red", lwd=3, lty=2)

plot(tajima$A.coluzzii_Mali_TD~tajima$stop,type="l",col="dodgerblue4",xlab = "Position", ylab = "Tajima's D",cex.axis=1.6,cex.lab = 1.3,main = "A.coluzzii(Mali)")
abline(h =c(0,qpre_td_ColM) ,col=c("black","dodgerblue4","sienna3","plum4"), lwd=1.5, lty=1:4)
#abline(v = c(0,24000000), col="red", lwd=3, lty=2)
#abline(v = c(start,end), col="red", lwd=3, lty=2)

plot(tajima$A.gambiae_Benin_TD~tajima$stop,type="l",col="sienna3",xlab = "Position", ylab = "Tajima's D",cex.axis=1.6,cex.lab = 1.3,main = "A.gambiae(Benin)")
abline(h =c(0,qpre2_td_GamB) ,col=c("black","sienna3"), lwd=1.5, lty=1:4)
#abline(v = c(0,24000000), col="red", lwd=3, lty=2)
#abline(v = c(start,end), col="red", lwd=3, lty=2)

plot(tajima$A.gambiae_Mali_TD~tajima$stop,type="l",col="plum4",xlab = "Position", ylab = "Tajima's D",cex.axis=1.6,cex.lab = 1.3,main = "A.gambiae(Mali)")
abline(h =c(0,qpre2_td_GamM) ,col=c("black","plum4"), lwd=1.5, lty=1:4)
#abline(v = c(0,24000000), col="red", lwd=3, lty=2)
#abline(v = c(start,end), col="red", lwd=3, lty=2)

dev.off()




#######################################################################
################### plot dxy ###################################
#######################################################################


setwd("~/Documents/malaria_vector_reserach/input_data/introgression_mali_benin/2L")
dxy<-read.table(file = "2L_sparrow_dxy_sites_all_pops_Benin_Mali.txt", header = TRUE)

the_quantile_td_p <- 0.99
dxy1 <- quantile (dxy[,"A.coluzzii_Benin.A.coluzzii_Mali_dxy"], the_quantile_td_p,na.rm = T)
dxy2<- quantile (dxy[,"A.coluzzii_Benin.A.gambiae_Benin_dxy"], the_quantile_td_p,na.rm = T)
dxy3 <- quantile (dxy[,"A.coluzzii_Benin.A.gambiae_Mali_dxy"], the_quantile_td,na.rm = T)
dxy4 <- quantile (dxy[,"A.coluzzii_Mali.A.gambiae_Benin_dxy"], the_quantile_td,na.rm = T)
dxy5 <- quantile (dxy[,"A.coluzzii_Mali.A.gambiae_Mali_dxy"], the_quantile_td,na.rm = T)
dxy6 <- quantile (dxy[,"A.gambiae_Benin.A.gambiae_Mali_dxy"], the_quantile_td,na.rm = T)

pdf(file = "2L_An.gambiae_An.gambiae_Dxy_sep.pdf", width = 14, height= 10)

par(mfrow=c(3,2))
plot(dxy$A.coluzzii_Benin.A.coluzzii_Mali_dxy~dxy$stop,type="l",col="darkcyan",xlab = "Position", cex.axis=1.6,cex.lab = 1.3,ylab = "dxy", main = "A.coluzzii(Benin)_A.coluzzii(Mali)")
abline(h =(dxy1) ,col=c("darkcyan"), lwd=0.7, lty=4)
abline(v = c(0,24000000), col="red", lwd=3, lty=2)
#abline(v = (2422651), col="red", lwd=3, lty=2)

plot(dxy$A.coluzzii_Benin.A.gambiae_Benin_dxy~dxy$stop,type="l",col="dodgerblue4",xlab = "Position", ylab = "dxy",cex.axis=1.6,cex.lab = 1.3,main = "A.coluzzii(Benin)_A.gambiae(Benin)")
abline(h =(dxy2) ,col=c("dodgerblue4"), lwd=0.7, lty=4)
abline(v = c(0,24000000), col="red", lwd=3, lty=2)
#abline(v = (2422651), col="red", lwd=3, lty=2)

plot(dxy$A.coluzzii_Benin.A.gambiae_Mali_dxy~dxy$stop,type="l",col="sienna3",xlab = "Position", ylab = "dxy",cex.axis=1.6,cex.lab = 1.3,main = "A.coluzzii(Benin)_A.gambiae(Mali)")
abline(h =(dxy3) ,col=c("sienna3"), lwd=0.7, lty=4)
abline(v = c(0,24000000), col="red", lwd=3, lty=2)
#abline(v = (2422651), col="red", lwd=3, lty=2)

plot(dxy$A.coluzzii_Mali.A.gambiae_Mali_dxy~dxy$stop,type="l",col="plum4",xlab = "Position", ylab = "Tajima's D",cex.axis=1.6,cex.lab = 1.3,main = "A.coluzzii(Mali)_A.gambiae(Mali)")
abline(h =(dxy4) ,col=c("plum4"), lwd=0.7, lty=4)
abline(v = c(0,24000000), col="red", lwd=3, lty=2)
#abline(v = (2422651), col="red", lwd=3, lty=2)


plot(dxy$A.coluzzii_Mali.A.gambiae_Benin_dxy~dxy$stop,type="l",col="chartreuse4",xlab = "Position", ylab = "dxy",cex.axis=1.6,cex.lab = 1.3,main = "A.coluzzii(Mali)_A.gambiae(Benin)")
abline(h =(dxy5) ,col=c("chartreuse4"), lwd=0.7, lty=4)
abline(v = c(0,24000000), col="red", lwd=3, lty=2)
#abline(v = (2422651), col="red", lwd=3, lty=2)

plot(dxy$A.gambiae_Benin.A.gambiae_Mali_dxy~dxy$stop,type="l",col="brown4",xlab = "Position", ylab = "dxy",cex.axis=1.6,cex.lab = 1.3,main = "A.gambiae(Benin)_A.gambiae(Mali)")
abline(h =(dxy4) ,col=c("brown4"), lwd=0.7, lty=4)
abline(v = c(0,24000000), col="red", lwd=3, lty=2)
#abline(v = (2422651), col="red", lwd=3, lty=2)

dev.off()