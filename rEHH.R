## tutorial: https://cran.r-project.org/web/packages/rehh/vignettes/rehh.html
#rm(list = ls())

library(rehh)
library(vcfR)  ## to read vcf file format
library(ape)



setwd("~/Documents/malaria_vector_reserach/input_data/introgression_mali_benin/data_March_2020/ehh_analysis")
hh<-data2haplohh(hap_file = "flower.phased.vcf.ehh.input.with.locus.id.gz",
                 #map_file = "2L.marker_information.inp",
                 polarize_vcf = FALSE)

#cat(readLines("2L.marker_information.inp"), sep = "\n")

res <- calc_ehh(hh, 
                mrk = "2L:2422652", 
                include_nhaplo = TRUE)

res$ihh
pdf(file = "An.gambiae_Benin_EHH1_KDR_updated.pdf", height= 8, width = 8)
par(mar = c(5, 3, 4.1, 3))
plot(res, lwd=3, xlim=c(2400000,2440000),cex.axis=1.6,cex.lab = 1.3,col=c("orangered3","green3")) 
dev.off()    



##3.3 The function calc_ehhs()

res_ehhs <- calc_ehhs(hh, 
                 mrk = "2L:2422652", 
                 include_nhaplo = TRUE)

pdf(file = "An.gambaie_Benin_EHH2_KDR_updated.pdf", height= 8, width = 8)
plot(res_ehhs,lwd=3,xlim=c(2400000,2430000),cex.axis=1.6,cex.lab = 1.3,col = "seagreen")
dev.off()


## bifurcation plot
furcation <- calc_furcation(hh,
                            mrk = "2L:2422652")

pdf(file = "An.gambiae_Benin_haplotype_furcation_KDR_updated.pdf", height= 8, width = 8)
par(mar = c(5, 3, 4.1, 3))
plot(furcation,lwd = 0.8, legend=NULL,xlim=c(2400000,2450000),cex.axis=1.6,cex.lab = 1.3)
dev.off()



pdf(file = "An.gambiae_Benin_haplotype_furcation2_ind_KDR_updated.pdf", height= 8, width = 8)
par(mar = c(4, 7, 1, 7))
plot(furcation,lwd = 0.9,hap.names = hap.names(hh), cex= 0.5,legend=NULL,xlim=c(2400000,2450000))
dev.off()




#Visualization of the haplotype structure
#The function plot.haplohh()
## sequences in rows and markers in columns.
hh_subset <- subset(hh, select.mrk = 20000:23500)

pdf(file = "An.gambiae_coluzzi_haplotype_matrix_network.pdf", height= 8, width = 8)
oldpar <- par(mar = c(3, 2, 2, 2) + 0.1)
plot(
  hh_subset,
  mrk = "2L:2422652",
  group_by_allele = TRUE,
  ignore.distance = TRUE,
  col = c(NA, "red"),
  linecol = c("lightblue", "lightpink"),
  mrk.col = "black",
  cex = 0.1,
  lwd=2,
  pos.lab.hap = "none",
  pos.lab.mrk = "none"
)
par(oldpar)
dev.off()

#### phylogenetics of 
newick <- as.newick(furcation,
                    allele = 1,
                    side = "left",
                    hap.names = hap.names(hh))

library(ape)
tree <- read.tree(text = newick)

pdf(file = "An.gambiae.Benin_haplotype_tree.pdf", height= 8, width = 8)
plot(tree, 
     cex = 0.7, 
     direction = "leftwards", 
     edge.color = "blue",
     underscore = TRUE,
     no.margin = TRUE)
dev.off()


#The functions calc_haplen() and plot.haplen() ### plot haplotype length
haplen <- calc_haplen(furcation)
haplen$mrk.name
head(haplen$haplen)

pdf(file = "An.gambiae.MAli_haplotype_length_KDR_updated.pdf", height= 8, width = 8)
par(mar = c(5, 2, 4, 8))
plot(haplen,lwd= 3)
dev.off()


################################
###The function scan_hh() {#scanhh}

scan <- scan_hh(hh)

# perform scan applying calc_ehh and calc_ehhs to each marker
slow_scan_hh <- function(haplohh) {
  # create empty vectors of size nmrk
  IHH_A <- IHH_D <- IES <- INES <- vector("numeric", nmrk(haplohh))
  # invoke calc_ehh and calc_ehhs for each marker
  for (i in 1:nmrk(haplohh)) {
    res <- calc_ehh(haplohh, mrk = i)
    IHH_A[i] <- res$ihh["IHH_A"]
    IHH_D[i] <- res$ihh["IHH_D"]
    res <- calc_ehhs(haplohh, mrk = i)
    IES[i] <- res$IES
    INES[i] <- res$INES
  }
  # create data frame (the return value of this function)
  data.frame(IHH_A = IHH_A, 
             IHH_D = IHH_D,
             IES = IES,
             INES = INES)
}
system.time(slow_scan <- slow_scan_hh(hh))



wgscan.ihs.cgu <- ihh2ihs(scan)
head(wgscan.ihs.cgu$frequency.class)

#Distribution of standardized values: the function distribplot()

### manhattan plot
manhattanplot(wgscan.ihs.cgu,
              pval = TRUE,
              threshold = 4,
              main = "p-value of iHS (CGU cattle breed)")


#Genome wide score plots: the function manhattan() of package qqman

# extract data frame from result list
ihs <- wgscan.ihs.cgu$ihs
# create new data frame
wgscan.cgu.ihs.qqman <- data.frame(
  CHR = as.integer(factor(ihs$CHR, 
                          levels = unique(ihs$CHR))),
  # chromosomes as integers
  BP = ihs$POSITION,         # base pairs
  P = 10**(-ihs$LOGPVALUE),  # transform back to p-values
  SNP = row.names(ihs)       # SNP names
)

library(qqman)
manhattan(wgscan.cgu.ihs.qqman,
          col = c("red","green"),
          chrlabs = unique(ihs$CHR),
          suggestiveline = 4,
          highlight = "F1205400",
          annotatePval = 0.0001, ylim = 10)


###6.2 The functions calc_haplen() and plot.haplen()
haplen <- calc_haplen(furcation)
haplen$mrk.name
haplen$position
haplen$xlim
pdf(file = "ac.Benin_haplotype_length_1575Y.pdf", height= 8, width = 8)
#plot(haplen,lwd=6, legend= NULL)
plot(haplen,lwd=6)
dev.off()



Genome wide score plots: the function manhattan() of package qqman

