#
rm(list = ls())
# compress and index the VCF
#bgzip snp.vcf
#tabix -p vcf snp.vcf.gz


## tutorial and manual
#https://cran.r-project.org/web/packages/PopGenome/PopGenome.pdf
#https://cran.r-project.org/web/packages/PopGenome/vignettes/Whole_genome_analyses_using_VCF_files.pdf
#https://markravinet.github.io/Chapter8.html


### annuus.annuus-argophyllus.chrom9
#seqlength <-c(chr9=200254352)
#bins<-tileGenome(seqlength,tilewidth=5000000,cut.last.tile.in.chrom = T)
#write.table(bins,file="chrom9_bins",col.names = FALSE, row.names = TRUE)


library(tidyverse)
library(PopGenome)
setwd("/share/lanzarolab/users/shaghayegh/data/malli_vs_benin/data_March2020/2L_sliced_vcfs_hom_kdr_sfs_phylo/1-3MB_phased")
sparrows <- readVCF("2L_kdr_hom_1_3MB.phased.vcf.gz",numcols=10000, tid="2L",
                    from=1, to= 2999964, parallel=FALSE, gffpath=FALSE)


get.sum.data(sparrows)
show.slots(sparrows)


# check for population data
sparrows@populations

sparrow_info <- read_delim("./sparrow_inds.txt", delim = "\t")
#ind	pop
#ANN0801	ANN_G1
#ANN0802	ANN_G1
#ANN0803	ANN_G1
#ANN0804	ANN_G1


# now get the data for the popultions
populations <- split(sparrow_info$ind, sparrow_info$pop)
write.table(populations,"sparrow_inds_ordered_from_popgenome.txt")

# now set 
sparrows <- set.populations(sparrows, populations, diploid = T)

# set chromosome size
chr1 <- 2999964   ## chromosme size should match with the "to= ..) in the read vcf file

# set window size and window jump
window_size <- 50000
window_jump <- 10000



#use seq to find the start points of each window
window_start <- seq(from = 1, to = chr1, by = window_jump)
# add the size of the window to each start point 
window_stop <- window_start + window_size


# no windows start before the end of chromosome 8
sum(window_start > chr1)
# but some window stop positions do occur past the final point
sum(window_stop > chr1)


# remove windows from the start and stop vectors
window_start <- window_start[which(window_stop < chr1)]
window_stop <- window_stop[which(window_stop < chr1)]

##this highlights an important point, our final window actually falls short of the end of the chromosome. You can check this like so:

chr1 - window_stop[length(window_stop)]



# save as a data.frame
windows <- data.frame(start = window_start, stop = window_stop, 
                      mid = window_start + (window_stop-window_start)/2)                     

# make a sliding window dataset
sparrows_sw <- sliding.window.transform(sparrows, width = 50000, jump = 10000, type = 2)


#### Extracting data ####
### calculate diversity.statistics (nucleotide diversity)
sparrows_sw <- diversity.stats(sparrows_sw, pi = TRUE)

### calculate diversity statistics (Fst,dxy)
sparrows_sw <- F_ST.stats(sparrows_sw, mode = "nucleotide")


#### calculate neutrality.statistics (tajima's_D, segregating sites)
sparrows_sw <- neutrality.stats(sparrows_sw)



# convert extrcated stats for the visualization
## extrcat nucleotide.diversity (pi) stats and correct for window size
nd <- sparrows_sw@nuc.diversity.within/50000    ## based on the window_size we chose
pops <- c("06BANA0012_gM","06BANA0013_gM","06BANA0016_gM","06DONE0104_gM","06FOUN0006_gM","06FOUN0007_gM","06FOUN0026_gM","06SAMA0130_gM","06SAMA0211_gM","06SELI0115_gM","06SELI0123_gM","12SELI0002_cM","12SELI0003_cM","12SELI0012_cM","12SELI0021_cM","12SELI0026_cM","2014ABO001_gB","2014ABO002_gB","2014ABO003_cB","2014ABO004_gB","2014ABO005_gB","2014ABO006_gB","2014COV001_cB","2014COV002_cB","2014COV003_cB","2014COV004_cB","2014COV006_cB","2014COV009_cB","2014COV012_cB","2014DONE0003_cM","2014DONE0032_cM","2014DONE0041_cM","2014DONE0050_cM","2014TOR001_gB","2014TOR002_cB","2014TOR004_cB","2014TOR005_gB" ,"2014TOR006_gB", "2014TOR007_cB" ,"2014TOR008_cB","2014TOR009_gB" ,"2014TOR010_cB","2014TOR011_cB","2015DUNA0099_cM","2015DUNA0100_cM","2015FANZ0013_cM","2015FANZ0043_cM","2015FANZ0044_cM")

colnames(nd) <- paste0(pops, "_pi")
sparrow_nd<-cbind(windows,nd)
sparrow_nd$mid<-(sparrow_nd[,3]/10^6)
write.table(sparrow_nd, file = "sparrow_nd_all_Inds_50k_10kstepping.txt", col.names = TRUE, row.names = FALSE,quote = FALSE, sep = "\t")



### Tajima's D and number of segregating sites
tajima <- sparrows_sw@Tajima.D
colnames(tajima) <- paste0(pops, "_TD")
sparrow_tajima<-cbind(windows,tajima)
sparrow_tajima$mid<-(sparrow_tajima[,3]/10^6)
write.table(sparrow_tajima, file = "sparrow_Tajimas_D_all_Inds_5k_1kstepping.txt", col.names = TRUE, row.names = FALSE,quote = FALSE, sep = "\t")


sg_sites<-sparrows_sw@n.segregating.sites
colnames(sg_sites) <- paste0(pops, "_Nseg")
sparrow_sg_sites<-cbind(windows,sg_sites)
sparrow_sg_sites$mid<-(sparrow_sg_sites[,3]/10^6)
write.table(sparrow_sg_sites, file = "sparrow_sparrow_sg_sites_all_INDS_5k_1kstepping.txt", col.names = TRUE, row.names = FALSE,quote = FALSE, sep = "\t")


### DXY
# extract dxy - pairwise absolute nucleotide diversity
dxy <- get.diversity(sparrows_sw, between = T)[[2]]/50000

### Fst
fst <- t(sparrows_sw@nuc.F_ST.pairwise)


# get column names right
x <- colnames(fst)
x <- colnames(dxy)
# replace all occurrences of pop1 with house
x <- sub("\\<pop1\\>", "06BANA0012_gM", x)
x <- sub("\\<pop2\\>", "06BANA0013_gM", x)
x <- sub("\\<pop3\\>", "06BANA0016_gM", x)
x <- sub("\\<pop4\\>", "06DONE0104_gM", x)
x <- sub("\\<pop5\\>", "06FOUN0006_gM", x)
x <- sub("\\<pop6\\>", "06FOUN0007_gM", x)
x <- sub("\\<pop7\\>", "06FOUN0026_gM", x)
x <- sub("\\<pop8\\>", "06SAMA0130_gM", x)
x <- sub("\\<pop9\\>", "06SAMA0211_gM", x)
x <- sub("\\<pop10\\>", "06SELI0115_gM", x)
x <- sub("\\<pop11\\>", "06SELI0123_gM", x)
x <- sub("\\<pop12\\>", "12SELI0002_cM", x)
x <- sub("\\<pop13\\>", "12SELI0003_cM", x)
x <- sub("\\<pop14\\>", "12SELI0012_cM", x)
x <- sub("\\<pop15\\>", "12SELI0021_cM", x)
x <- sub("\\<pop16\\>", "12SELI0026_cM", x)
x <- sub("\\<pop17\\>", "2014ABO001_gB", x)
x <- sub("\\<pop18\\>", "2014ABO002_gB", x)
x <- sub("\\<pop19\\>", "2014ABO003_cB", x)
x <- sub("\\<pop20\\>", "2014ABO004_gB", x)
x <- sub("\\<pop21\\>", "2014ABO005_gB", x)
x <- sub("\\<pop22\\>", "2014ABO006_gB", x)
x <- sub("\\<pop23\\>", "2014COV001_cB", x)
x <- sub("\\<pop24\\>", "2014COV002_cB", x)
x <- sub("\\<pop25\\>", "2014COV003_cB", x)
x <- sub("\\<pop26\\>", "2014COV004_cB", x)
x <- sub("\\<pop27\\>", "2014COV006_cB", x)
x <- sub("\\<pop28\\>", "2014COV009_cB", x)
x <- sub("\\<pop29\\>", "2014COV012_cB", x)
x <- sub("\\<pop30\\>", "2014DONE0003_cM", x)
x <- sub("\\<pop31\\>", "2014DONE0032_cM", x)
x <- sub("\\<pop32\\>", "2014DONE0041_cM", x)
x <- sub("\\<pop33\\>", "2014DONE0050_cM", x)
x <- sub("\\<pop34\\>", "2014TOR001_gB", x)
x <- sub("\\<pop35\\>", "2014TOR002_cB", x)
x <- sub("\\<pop36\\>", "2014TOR004_cB", x)
x <- sub("\\<pop37\\>", "2014TOR005_gB", x)
x <- sub("\\<pop38\\>", "2014TOR006_gB", x)
x <- sub("\\<pop39\\>", "2014TOR007_cB", x)
x <- sub("\\<pop40\\>", "2014TOR008_cB", x)
x <- sub("\\<pop41\\>", "2014TOR009_gB", x)
x <- sub("\\<pop42\\>", "2014TOR010_cB", x)
x <- sub("\\<pop43\\>", "2014TOR011_cB", x)
x <- sub("\\<pop44\\>", "2015DUNA0099_cM", x)
x <- sub("\\<pop45\\>", "2015DUNA0100_cM", x)
x <- sub("\\<pop46\\>", "2015FANZ0013_cM", x)
x <- sub("\\<pop47\\>", "2015FANZ0043_cM", x)
x <- sub("\\<pop48\\>", "2015FANZ0044_cM", x)


# replace forward slash
x <- sub("/", ":", x)
# look at x to confirm the replacement has occurred



#paste0(x, "_fst")
#paste0(x, "_dxy")

colnames(fst) <- paste0(x, "")
sparrow_fst<-cbind(windows,fst)
sparrow_fst$mid<-(sparrow_fst[,3]/10^6)


colnames(dxy) <- paste0(x, "")
sparrow_dxy<-cbind(windows,dxy)
sparrow_dxy$mid<-(sparrow_dxy[,3]/10^6)


write.table(sparrow_fst, file = "sparrow_sparrow_fst_all_INDS48_5k_1kstepping.txt", col.names = TRUE, row.names = FALSE,quote = FALSE, sep = "\t")
write.table(sparrow_dxy, file = "sparrow_sparrow_dxy_all_INDS48_5k_1kstepping.txt", col.names = TRUE, row.names = FALSE,quote = FALSE, sep = "\t")


dxy_mean<-data.frame(colMeans(sparrow_dxy))

dxy_mean <- data.frame(names = row.names(dxy_mean), dxy_mean)
rownames(dxy_mean)<-NULL
dxy_mean<-dxy_mean[-c(1,2,3),]
colnames(dxy_mean)<-c("comp","dxy")


#### make a file to include other way comparison
other<-dxy_mean[,1]
other_c1<- gsub(":.*$", "", dxy_mean$comp)
other_c2 <- sub('.*:\\s*', '', dxy_mean$comp)
other<-data.frame(cbind(other_c1,other_c2))
other$comp<-paste(other$other_c2,other$other_c1, sep = ":")
other$comp_present<-paste(other$other_c1,other$other_c2, sep = ":")

other_ad<-merge(other,dxy_mean,by.x = "comp_present",by.y = "comp")
other_ad<-other_ad[,c(4:5)]

### load ind file with 0 dxy within each individual
ind<-read.table(file="sparrow_inds_ordered_from_popgenome.txt", header = FALSE)
#06BANA0012_gM	0
#06BANA0013_gM	0
#06BANA0016_gM	0
#06DONE0104_gM	0
#06FOUN0006_gM	0

ind$V3<-paste(ind$V1,ind$V1, sep = ":")
ind<-ind[,c(3,2)]
colnames(ind)<-c("comp","dxy")

dxy_merged <- do.call("rbind", list(dxy_mean,other_ad,ind))
dxy_merged$ind1<- gsub(":.*$", "", dxy_merged$comp)
dxy_merged$ind2<- gsub('.*:\\s*', "", dxy_merged$comp)
dxy_merged_format<-dxy_merged[,c(3,4,2)]

#ind1          ind2          dxy
#1 06BANA0012_gM 06BANA0013_gM 0.0005660339
#2 06BANA0012_gM 06BANA0016_gM 0.0007592203
#3 06BANA0012_gM 06DONE0104_gM 0.0006247458

## write the completed table of both-ways comparisons
write.table(dxy_merged_format,file = "sparrow_sparrow_dxy_all_INDS48_5k_1kstepping_adjusted_both_comparisons.txt",col.names = TRUE, row.names = FALSE,quote = FALSE, sep = "\t")

###
# setwd("~/Documents/malaria_vector_reserach/input_data/introgression_mali_benin/data_March_2020/1-3MB_kdr_hom/1-3MB_2L_phased_50K_10Kstepping_inds")
dxy_merged_format<-read.table(file = "sparrow_sparrow_dxy_all_INDS48_5k_1kstepping_adjusted_both_comparisons.txt", header = TRUE)

## create qunique charcater type ID
ch1<- gsub('.*_\\s*', "", dxy_merged_format$ind1)
ch2<-gsub('.*_\\s*', "", dxy_merged_format$ind2)
charcaters<-data.frame(cbind(ch1,ch2))
charcaters$unique<-paste(charcaters$ch1,charcaters$ch2,sep="_")


### open ind names with desired order
ind_des<-read.table(file ="individual_list_desired.txt", header = FALSE)
dxy_merged_format[order(ind_des$V1),]

dxy_merged_format<-cbind(dxy_merged_format,charcaters$unique)

comparisons<-c("gM_gM","gM_gB","gM_cM","gM_cB","gB_gM","gB_gB","gB_cM","gB_cB","cM_gM","cM_gB","cM_cM","cM_cB","cB_gM","cB_gB","cB_cM","cB_cB")

gM_gM<-dxy_merged_format[dxy_merged_format$`charcaters$unique`=="gM_gM",]
gM_gB<-dxy_merged_format[dxy_merged_format$`charcaters$unique`=="gM_gB",]
gM_cM<-dxy_merged_format[dxy_merged_format$`charcaters$unique`=="gM_cM",]
gM_cB<-dxy_merged_format[dxy_merged_format$`charcaters$unique`=="gM_cB",]
g1<-do.call("rbind", list(gM_gM,gM_gB,gM_cM,gM_cB))

gB_gM<-dxy_merged_format[dxy_merged_format$`charcaters$unique`=="gB_gM",]
gB_gB<-dxy_merged_format[dxy_merged_format$`charcaters$unique`=="gB_gB",]
gB_cM<-dxy_merged_format[dxy_merged_format$`charcaters$unique`=="gB_cM",]
gB_cB<-dxy_merged_format[dxy_merged_format$`charcaters$unique`=="gB_cB",]
g2<-do.call("rbind", list(gB_gM,gB_gB,gB_cM,gB_cB))


cM_gM<-dxy_merged_format[dxy_merged_format$`charcaters$unique`=="cM_gM",]
cM_gB<-dxy_merged_format[dxy_merged_format$`charcaters$unique`=="cM_gB",]
cM_cB<-dxy_merged_format[dxy_merged_format$`charcaters$unique`=="cM_cB",]
cM_cM<-dxy_merged_format[dxy_merged_format$`charcaters$unique`=="cM_cM",]
g3<-do.call("rbind", list(cM_gM,cM_gB,cM_cB,cM_cM))


cB_gM<-dxy_merged_format[dxy_merged_format$`charcaters$unique`=="cB_gM",]
cB_gB<-dxy_merged_format[dxy_merged_format$`charcaters$unique`=="cB_gB",]
cB_cM<-dxy_merged_format[dxy_merged_format$`charcaters$unique`=="cB_cM",]
cB_cB<-dxy_merged_format[dxy_merged_format$`charcaters$unique`=="cB_cB",]
g4<-do.call("rbind", list(cB_gM,cB_gB,cB_cM,cB_cB))


dxy_merged_format<- do.call("rbind", list(g1,g2,g3,g4))

## desired order of individuals and reorder the main sxy table based on the desired table
#desired<-read.table(file = "individual_list_desired.txt", header = FALSE)
#dxy_merged_format<-dxy_merged_format[order(match(dxy_merged_format[,1],desired[,1])),]

dxy_merged_format<-dxy_merged_format[,c(1:3)]

## make pivot table ##
dxy_merged_format$dxy<-as.numeric(as.character(dxy_merged_format$dxy))
sata_pivot_r<-reshape2::dcast(dxy_merged_format, ind1 ~ ind2, value.var="dxy", fun.aggregate=sum)


rnames_r <- sata_pivot_r[,1]
mat_data_r <- data.matrix(sata_pivot_r[,2:ncol(sata_pivot_r)])
rownames(mat_data_r) <- rnames_r

heatmap.2(mat_data_r,trace="none",cexRow=0.3,cexCol =0.5,srtCol=45,cex.main=0.01,margins = c(8,8))



#mat_data_r<-mat_data_r[order(mat_data_r[,ncol(mat_data_r)],decreasing=T),]
my_palette <- colorRampPalette(c("red", "pink", "blue"))(n = 299)
#col_breaks = c(seq(-4,-1,length=100),  # for red
#               seq(-1,1,length=100),              # for yellow
#               seq(1,4,length=100))    

pdf(file = "2L_1-3MB_dxy_heatmap.pdf", height= 12, width=12)
heatmap.2(mat_data_r,col=my_palette,trace="none",margins = c(8,8))
dev.off()

my_palette <- colorRampPalette(c("red", "pink", "blue"))(n = 229)
heatmap.2(mat_data_r,col=my_palette,trace="none",margins = c(8,8))

