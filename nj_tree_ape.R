
## nj tree analysis in R (ape package)
library(adegenet)
library(ape)
## tutorial: https://fuzzyatelin.github.io/bioanth-stats/module-24/module-24.html
## introduction to phylogenetics in r: http://www.phytools.org/Cordoba2017/ex/2/Intro-to-phylogenies.html
setwd("~/Documents/malaria_vector_reserach/input_data/introgression_mali_benin/data_March_2020/2-3MB_kdr_hom/phylogenetics/phylo_input_from_IR_perl")

dna <- read.FASTA(file= "clustalo-I20200409-192236-0447-62229962-p1m_fasta", type = "DNA")   ### vcf converted into FASTA

## load annotation file 
my_annot <- read.table(file = "fasta_sample_ordered_annotaion.txt", header=FALSE)
colnames(my_annot)<-c("sample_ID","cluster")
#sample_ID cluster
#1 06FOUN0007   gam_M
#2 12SELI0012   col_M


#Step 1: find genetic distances for pairs of individuals (in our case, isolates)
D <- dist.dna(dna, model = "TN93")
D <- dist.dna(dna)
length(D) #number of pairwise distances, computed as n(n-1)/2

## plot matrix
temp <- as.data.frame(as.matrix(D))
table.paint(temp, cleg=0, clabel.row=.5, clabel.col=.5) #darker shades of gray mean a larger dist


temp <- t(as.matrix(D))
temp <- temp[,ncol(temp):1]


## step 2: make nj tree
tre <- nj(D)
class(tre)
## [1] "phylo"
tre <- ladderize(tre)
tre

## Unrooted; includes branch lengths.
plot(tre, cex=.7)
title("NJ tree; 2-3 MB (2L)")

## plot with color codes
plot(tre, cex=.7,show.tip=FALSE)

title("Unrooted NJ tree")
myPal <- colorRampPalette(c("red","darkblue","darkgreen"))
tiplabels(my_annot$cluster,bg=num2col(my_annot$year, col.pal=myPal),adj = c(0, 0),cex=.7) #we use the annot dataset to get our years
#temp <- pretty(1993:2008, 5)
#legend("bottomleft", fill=num2col(temp, col.pal=myPal), leg=temp, ncol=2)


# or 
h_cluster <- hclust(D, method = "average", members = NULL) # method = average is used for UPGMA, members can be equal to NULL or a vector with a length of size D
plot(h_cluster, cex = 0.6)
#title("NJ tree; 2-3 MB (2L)")



plot(tre, show.tip=FALSE) # gets rid of the labels on the end, refer to the first tree depicted above
title("Unrooted NJ tree")
#myPal <- colorRampPalette(c("red","yellow","green","blue"))
#tiplabels(annot$species, bg=num2col(annot$species, col.pal=myPal), cex=.5) #we use the annot dataset to get our years
#temp <- pretty(1993:2008, 5)
#legend("bottomleft", fill=num2col(temp, col.pal=myPal), leg=temp, ncol=2)


plot(tre, type="unrooted", show.tip=FALSE)
title("Unrooted NJ tree")
tiplabels(tre$tip.label, 
          cex=.5)


## cluster dandrogram

h_cluster <- hclust(D, method = "average", members = NULL) # method = average is used for UPGMA, members can be equal to NULL or a vector with a length of size D
plot(h_cluster, cex = 0.6)


###
##upgama
tre3 <- as.phylo(hclust(D,method="average"))
y <- as.vector(as.dist(cophenetic(tre3)))
plot(x, y, xlab="original pairwise distances", ylab="pairwise distances on the tree",
     main="Is UPGMA appropriate?", pch=20, col=transp("black",.1), cex=3)
abline(lm(y~x), col="red")


cor(x,y)^2

plot(tre3)


myBoots <- boot.phylo(tre2, dna, function(e)
  root(nj(dist.dna(e, model = "TN93")),1))
myBoots


########
##Bootstrapping

myBoots <- boot.phylo(tre2, dna, function(e) root(nj(dist.dna(e, model = "TN93")),1))







#######



