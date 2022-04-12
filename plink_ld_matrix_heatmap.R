
#!/usr/bin/env Rscript
### this script use LDheatmap to plots the plink ld matrix 

#rm(list = ls())

library(LDheatmap)
library(snpStats)
library(Matrix)
setwd("/share/lanzarolab/users/shaghayegh/data/malli_vs_benin/data_March2020/2L_sliced_vcfs_hom_kdr/2-3MB_Unphased/ld_analysis_with_plink")

colB<-read.table(file = "col_Benin_2L_2-3MB_maf0.03..ld")
ld_colB<-as.matrix(colB)
colM<-read.table(file = "col_Mali_2L_2-3MB_maf0.03..ld")
ld_colM<-as.matrix(colM)
gamB<-read.table(file = "gam_Benin_2L_2-3MB_maf0.03..ld")
ld_gamB<-as.matrix(gamB)
gamM<-read.table(file = "gam_Mali_2L_2-3MB_maf0.03..ld")
ld_gamM<-as.matrix(gamM)



#mfrow=c(2, 2)
jpeg(file = "2L_kdr_hom_1_3MB_coluzzii_Benin.maf0.03.jpg", height = 800, width = 800)
color_spectrum = colorRampPalette(c("Red", "Yellow")) #Creat a color spectrum
LDheatmap(ld_colB,
          color = color_spectrum(5),add.map=FALSE, title = "A.colluzzii-Benin-late 2L:2-3Mbp")
dev.off()



jpeg(file = "2L_kdr_hom_1_3MB_coluzzii_Mali.maf0.03.jpg", height = 800, width = 800)
color_spectrum = colorRampPalette(c("Red", "Yellow")) #Creat a color spectrum
LDheatmap(ld_colM,
          color = color_spectrum(5),add.map=FALSE, title = "A.colluzzii-Malin-late 2L:2-3Mbp")
dev.off()


jpeg(file = "2L_kdr_hom_1_3MB_gambiae_Benin.maf0.03.jpg", height = 800, width = 800)
color_spectrum = colorRampPalette(c("Red", "Yellow")) #Creat a color spectrum
LDheatmap(ld_gamB,
          color = color_spectrum(5),add.map=FALSE, title = "A.gambiae-Benin 2L:2-3Mbp")
dev.off()


jpeg(file = "2L_kdr_hom_1_3MB_gambiae_Mali.maf0.03.jpg", height = 800, width = 800)
color_spectrum = colorRampPalette(c("Red", "Yellow")) #Creat a color spectrum
LDheatmap(ld_gamM,
          color = color_spectrum(5),add.map=FALSE, title = "A.gambiae-Mali 2L:2-3Mbp")
dev.off()






