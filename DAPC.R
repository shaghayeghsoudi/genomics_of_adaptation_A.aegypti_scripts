
rm(list = ls())
library(adegenet)
library(vcfR)
library(poppr)
library(ape)
library(RColorBrewer)

## DAPC
### tutorial: https://grunwaldlab.github.io/Population_Genetics_in_R/gbs_analysis.html
setwd("~/Documents/malaria_vector_reserach/input_data/introgression_mali_benin/data_March_2020/2-3MB_kdr_hom/DAPC_2-3MB")
rubi.VCF <- read.vcfR("2L_kdr_hom_2_3MB_phased.vcf")
rubi.VCF


pop.data <- read.table("population_data.gbs.txt", sep = "\t", header = TRUE)

#We can now check that all the samples in the VCF and the population data frame are included:
all(colnames(rubi.VCF@gt)[-1] == pop.data$AccessID)


#Converting the dataset to a genlight object
gl.rubi <- vcfR2genlight(rubi.VCF)
pop(gl.rubi) <- pop.data$Country


###Distance matrices and plot 
tree <- aboot(gl.rubi, tree = "upgma", distance = bitwise.dist, sample = 48, showtree = F, cutoff = 60, quiet = T)

pdf(file = "genetic_distance_random20K_48_samples_introgressed_exclud.pdf", width = 15, height= 15)
cols <- brewer.pal(n = nPop(gl.rubi), name = "Dark2")
plot.phylo(tree, cex = 1, font = 2, adj = 0, tip.color =  cols[pop(gl.rubi)])
nodelabels(tree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.6,font = 2, xpd = TRUE)
#legend(35,10,c("A.colluzzi(Benin)","A.colluzzi(Mali)","A.gambiae(Benin)","A.gambiae(Mali)")),cols, border = FALSE, bty = "n")
legend("topright", legend = c("A.colluzzi(Benin)","A.colluzzi(Mali)","A.gambiae(Benin)","A.gambiae(Mali)"), fill = cols, border = FALSE, bty = "n", cex = 0.9)

axis(side = 1,cex=2.5)
title(xlab = "Genetic distance (proportion of loci that are different)")
dev.off()


### plot Minimum spanning networks
library(igraph)

rubi.dist <- bitwise.dist(gl.rubi)
rubi.msn <- poppr.msn(gl.rubi, rubi.dist, showplot = FALSE, include.ties = T)

node.size <- rep(3, times = nInd(gl.rubi))
names(node.size) <- indNames(gl.rubi)
vertex.attributes(rubi.msn$graph)$size <- node.size

set.seed(9)
pdf(file = "genetic_distance_network_random20K_48_samples_introgressed_exclud.pdf", width = 10, height= 10)
plot_poppr_msn(gl.rubi, rubi.msn , palette = brewer.pal(n = nPop(gl.rubi), name = "Dark2"), gadj = 70)
dev.off()


#Principal components analysis
rubi.pca <- glPca(gl.rubi, nf = 3)
barplot(100*rubi.pca$eig/sum(rubi.pca$eig), col = heat.colors(50), main="PCA Eigenvalues")
title(ylab="Percent of variance\nexplained", line = 2)
title(xlab="Eigenvalues", line = 1)
rubi.pca.scores <- as.data.frame(rubi.pca$scores)
rubi.pca.scores$pop <- pop(gl.rubi)

library(ggplot2)
set.seed(9)
p <- ggplot(rubi.pca.scores, aes(x=PC1, y=PC2, colour=pop)) 
p <- p + geom_point(size=3)
p <- p + stat_ellipse(level = 0.95, size = 1)
p <- p + scale_color_manual(values = cols) 
p <- p + geom_hline(yintercept = 0) 
p <- p + geom_vline(xintercept = 0) 
p <- p + theme_bw()
p


pnw.dapc <- dapc(gl.rubi, n.pca = 4, n.da = 48)
scatter(pnw.dapc, col = cols, cex = 2, legend = TRUE, clabel = F, posi.leg = "bottomleft", scree.pca = TRUE,
        posi.pca = "topleft", cleg = 0.75,pch = 17:22)        


myCol <- c("darkblue","purple","orange","red")
scatter(pnw.dapc,  posi.da="bottomright", bg="white",
        pch=17:22, cstar=0, col=myCol, scree.pca=TRUE,cex= 3,
        posi.pca="bottomleft")


