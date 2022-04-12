
rm(list = ls())
setwd("~/Documents/malaria_vector_reserach/input_data/introgression_mali_benin/data_March_2020/pop.genome.out.30k.3k")

## analyze nucleotide diversity data
dat_2R<-read.table(file = "sparrow_nd_2R_48_inds_30K_3K.txt", header = TRUE)
dat_2R$chrom<-rep("2R")
dat_2R$status<-rep("non_intogressed")

dat_3L<-read.table(file = "sparrow_nd_3L_48_inds_30K_3K.txt", header = TRUE)
dat_3L$chrom<-rep("3L")
dat_3L$status<-rep("non_intogressed")

dat_3R<-read.table(file = "sparrow_nd_3R_48_inds_30K_3K.txt", header = TRUE)
dat_3R$chrom<-rep("3R")
dat_3R$status<-rep("non_intogressed")

all_three<- do.call("rbind", list(dat_2R,dat_3L,dat_3R))
#d<-density(all_three$col_B_pi)
#plot(d, xlim = c(-0.00001,0.0003))

## load 2L and limit to the intogressed 1-3BM region
dat_2L<-read.table(file = "sparrow_nd_2L_48_inds_30K_3K.txt", header = TRUE)
dat_2L$chrom<-rep("2L")
dat_2L<-dat_2L[dat_2L$stop<3000000,]
dat_2L$status<-rep("intogressed")


### limit the data to 2-3MB
#dat_2L_lim<-dat_2L[dat_2L$start>2000000 & dat_2L$stop<3000000,]

data_all<-do.call("rbind", list(all_three,dat_2L))
#data_all<-do.call("rbind", list(all_three,dat_2L_lim))
## density plot
data_all$status<-factor(data_all$status)

p1<-ggplot(data_all,aes(x=col_B_pi,fill=status))+geom_density(alpha=.8)+ xlim(-0.00001,0.0003)+ ylim(0,80000)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.text=element_text(size=15),
       axis.title=element_text(size=14))
ggsave(filename="whole_genome_vs_2L_2-3MB_pi_density_ColB.pdf", plot=p1,width = 8, height = 8)


p2<-ggplot(data_all,aes(x=col_M_pi,fill=status))+geom_density(alpha=.8)+ xlim(-0.00001,0.0003)+ ylim(0,80000)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
       panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.text=element_text(size=15),
       axis.title=element_text(size=14))
ggsave(filename="whole_genome_vs_2L_2-3MB_pi_density_ColM.pdf", plot=p2,width = 8, height = 8)



p3<-ggplot(data_all,aes(x=gam_B_pi,fill=status))+geom_density(alpha=.8)+ xlim(-0.00001,0.0003)+ ylim(0,80000)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.text=element_text(size=15),
      axis.title=element_text(size=14))
ggsave(filename="whole_genome_vs_2L_2-3MB_pi_density_GamB.pdf", plot=p3,width = 8, height = 8)



p4<-ggplot(data_all,aes(x=gam_M_pi,fill=status))+geom_density(alpha=.8)+ xlim(-0.00001,0.0003)+ ylim(0,80000)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
     panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.text=element_text(size=15),
     axis.title=element_text(size=14))
ggsave(filename="whole_genome_vs_2L_2-3MB_pi_density_GamM.pdf", plot=p4,width = 8, height = 8)


###########################################
###### analyze tajimas D ###
###########################################
dat_2R<-read.table(file = "sparrow_Tajimas_D_2R_48_inds_30K_3K.txt", header = TRUE)
dat_2R$chrom<-rep("2R")
dat_2R$status<-rep("non_intogressed")

dat_3L<-read.table(file = "sparrow_Tajimas_D_3L_48_inds_30K_3K.txt", header = TRUE)
dat_3L$chrom<-rep("3L")
dat_3L$status<-rep("non_intogressed")

dat_3R<-read.table(file = "sparrow_Tajimas_D_3R_48_inds_30K_3K.txt", header = TRUE)
dat_3R$chrom<-rep("3R")
dat_3R$status<-rep("non_intogressed")

all_three<- do.call("rbind", list(dat_2R,dat_3L,dat_3R))
#d<-density(all_three$col_B_pi)
#plot(d, xlim = c(-0.00001,0.0003))

## load 2L and limit to the intogressed 1-3BM region
dat_2L<-read.table(file = "sparrow_Tajimas_D_2L_48_inds_30K_3K.txt", header = TRUE)
dat_2L$chrom<-rep("2L")
dat_2L<-dat_2L[dat_2L$stop<3000000,]
dat_2L$status<-rep("intogressed")


### limit the data to 2-3MB
#dat_2L_lim<-dat_2L[dat_2L$start>2400000 & dat_2L$stop<2500000,]

data_all<-do.call("rbind", list(all_three,dat_2L))
#data_all<-do.call("rbind", list(all_three,dat_2L_lim))
## density plot
data_all$status<-factor(data_all$status)

p1<-ggplot(data_all,aes(x=col_B_TD,fill=status))+geom_density(alpha=.8)+ xlim(-3,+4)+ ylim(0,0.8)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                                     panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.text=element_text(size=15),
                                                                                                                                                                                                           axis.title=element_text(size=14))
ggsave(filename="whole_genome_vs_2L_2-3MB_TjD_density_ColB.pdf", plot=p1,width = 8, height = 8)



#colm<-data_all[,c(1,2,3,5,8,9)]
#colm<-na.omit(colm)
p2<-ggplot(data_all,aes(x=col_M_TD,fill=status))+geom_density(alpha=.8)+ xlim(-3,+4)+ ylim(0,0.8)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                                    panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.text=element_text(size=15),
                                                                                                                                                                                                          axis.title=element_text(size=14))
ggsave(filename="whole_genome_vs_2L_2-3MB_TjD_density_ColM.pdf", plot=p2,width = 8, height = 8)



p3<-ggplot(data_all,aes(x=gam_B_TD,fill=status))+geom_density(alpha=.8)+ xlim(-3,+4)+ ylim(0,0.8)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                                    panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.text=element_text(size=15),
                                                                                                                                                                                                          axis.title=element_text(size=14))
ggsave(filename="whole_genome_vs_2L_2-3MB_TjD_density_GamB.pdf", plot=p3,width = 8, height = 8)



p4<-ggplot(data_all,aes(x=gam_M_TD,fill=status))+geom_density(alpha=.8)+ xlim(-3,+4)+ ylim(0,0.8)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                                    panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.text=element_text(size=15),
                                                                                                                                                                                                          axis.title=element_text(size=14))
ggsave(filename="whole_genome_vs_2L_2-3MB_TjD_density_GamM.pdf", plot=p4,width = 8, height = 8)


############
### analyze FuLi-D


dat_2R<-read.table(file = "sparrow_FuliD_2R_48_inds_30K_3K.txt", header = TRUE)
dat_2R$chrom<-rep("2R")
dat_2R$status<-rep("non_intogressed")

dat_3L<-read.table(file = "sparrow_FuliD_3L_48_inds_30K_3K.txt", header = TRUE)
dat_3L$chrom<-rep("3L")
dat_3L$status<-rep("non_intogressed")

dat_3R<-read.table(file = "sparrow_FuliD_3R_48_inds_30K_3K.txt", header = TRUE)
dat_3R$chrom<-rep("3R")
dat_3R$status<-rep("non_intogressed")

all_three<- do.call("rbind", list(dat_2R,dat_3L,dat_3R))
#d<-density(all_three$col_B_pi)
#plot(d, xlim = c(-0.00001,0.0003))

## load 2L and limit to the intogressed 1-3BM region
dat_2L<-read.table(file = "sparrow_FuliD_2L_48_inds_30K_3K.txt", header = TRUE)
dat_2L$chrom<-rep("2L")
dat_2L<-dat_2L[dat_2L$stop<3000000,]
dat_2L$status<-rep("intogressed")


### limit the data to 2-3MB
#dat_2L_lim<-dat_2L[dat_2L$start>2000000 & dat_2L$stop<3000000,]

data_all<-do.call("rbind", list(all_three,dat_2L))
#data_all<-do.call("rbind", list(all_three,dat_2L_lim))
## density plot
data_all$status<-factor(data_all$status)

p1<-ggplot(data_all,aes(x=col_B_FuliD,fill=status))+geom_density(alpha=.8)+ xlim(-4,+4)+ ylim(0,0.8)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                         panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.text=element_text(size=15),
                                                                                                                                                                                               axis.title=element_text(size=14))
ggsave(filename="whole_genome_vs_2L_1-3MB_FuliD_density_ColB.pdf", plot=p1,width = 8, height = 8)



#colm<-data_all[,c(1,2,3,5,8,9)]
#colm<-na.omit(colm)
p2<-ggplot(data_all,aes(x=col_M_FuliD,fill=status))+geom_density(alpha=.8)+ xlim(-4,+4)+ ylim(0,0.8)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                        panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.text=element_text(size=15),
                                                                                                                                                                                              axis.title=element_text(size=14))
ggsave(filename="whole_genome_vs_2L_1-3MB_FuliD_density_ColM.pdf", plot=p2,width = 8, height = 8)



p3<-ggplot(data_all,aes(x=gam_B_FuliD,fill=status))+geom_density(alpha=.8)+ xlim(-4,+4)+ ylim(0,0.8)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                        panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.text=element_text(size=15),
                                                                                                                                                                                              axis.title=element_text(size=14))
ggsave(filename="whole_genome_vs_2L_1-3MB_FuliD_density_GamB.pdf", plot=p3,width = 8, height = 8)



p4<-ggplot(data_all,aes(x=gam_M_FuliD,fill=status))+geom_density(alpha=.8)+ xlim(-4,+4)+ ylim(0,0.8)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                        panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.text=element_text(size=15),
                                                                                                                                                                                              axis.title=element_text(size=14))
ggsave(filename="whole_genome_vs_2L_1-3MB_FuliD_density_GamM.pdf", plot=p4,width = 8, height = 8)


##################################
########## analyze FuLi-F ######## 


dat_2R<-read.table(file = "sparrow_sparrow_FuliF_2R_48_inds_30K_3K.txt", header = TRUE)
dat_2R$chrom<-rep("2R")
dat_2R$status<-rep("non_intogressed")

dat_3L<-read.table(file = "sparrow_sparrow_FuliF_3L_48_inds_30K_3K.txt", header = TRUE)
dat_3L$chrom<-rep("3L")
dat_3L$status<-rep("non_intogressed")

dat_3R<-read.table(file = "sparrow_sparrow_FuliF_3R_48_inds_30K_3K.txt", header = TRUE)
dat_3R$chrom<-rep("3R")
dat_3R$status<-rep("non_intogressed")

all_three<- do.call("rbind", list(dat_2R,dat_3L,dat_3R))
#d<-density(all_three$col_B_pi)
#plot(d, xlim = c(-0.00001,0.0003))

## load 2L and limit to the intogressed 1-3BM region
dat_2L<-read.table(file = "sparrow_sparrow_FuliF_2L_48_inds_30K_3K.txt", header = TRUE)
dat_2L$chrom<-rep("2L")
dat_2L<-dat_2L[dat_2L$stop<3000000,]
dat_2L$status<-rep("intogressed")


### limit the data to 2-3MB
#dat_2L_lim<-dat_2L[dat_2L$start>2000000 & dat_2L$stop<3000000,]

data_all<-do.call("rbind", list(all_three,dat_2L))
#data_all<-do.call("rbind", list(all_three,dat_2L_lim))
## density plot
data_all$status<-factor(data_all$status)

p1<-ggplot(data_all,aes(x=col_B_Fulif,fill=status))+geom_density(alpha=.8)+ xlim(-4,+4)+ ylim(0,0.8)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                            panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.text=element_text(size=15),
                                                                                                                                                                                                  axis.title=element_text(size=14))
ggsave(filename="whole_genome_vs_2L_2-3MB_FuliF_density_ColB.pdf", plot=p1,width = 8, height = 8)



#colm<-data_all[,c(1,2,3,5,8,9)]
#colm<-na.omit(colm)
p2<-ggplot(data_all,aes(x=col_M_Fulif,fill=status))+geom_density(alpha=.8)+ xlim(-4,+4)+ ylim(0,0.8)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                           panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.text=element_text(size=15),
                                                                                                                                                                                                 axis.title=element_text(size=14))
ggsave(filename="whole_genome_vs_2L_2-3MB_FuliF_density_ColM.pdf", plot=p2,width = 8, height = 8)



p3<-ggplot(data_all,aes(x=gam_B_Fulif,fill=status))+geom_density(alpha=.8)+ xlim(-4,+4)+ ylim(0,0.8)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                           panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.text=element_text(size=15),
                                                                                                                                                                                                 axis.title=element_text(size=14))
ggsave(filename="whole_genome_vs_2L_2-3MB_FuliF_density_GamB.pdf", plot=p3,width = 8, height = 8)



p4<-ggplot(data_all,aes(x=gam_M_Fulif,fill=status))+geom_density(alpha=.8)+ xlim(-4,+4)+ ylim(0,0.8)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                           panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.text=element_text(size=15),
                                                                                                                                                                                                 axis.title=element_text(size=14))
ggsave(filename="whole_genome_vs_2L_2-3MB_FuliF_density_GamM.pdf", plot=p4,width = 8, height = 8)


