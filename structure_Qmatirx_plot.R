
setwd("~/Documents/malaria_vector_reserach/input_data/structure/fastStructure")
library(pophelper)
library(gridExtra)
library(mapplots)

# convert q-matrix run files to R qlist object
sfiles <- list.files("output_file_entire_ld.pruned_RUN1", pattern="*.meanQ", full.names = TRUE)
slist <- readQ(files=sfiles)

#readQ(files=sfiles,filetype="structure")
attributes(slist)
attributes(slist[[1]])

slist <- readQ(files=sfiles, indlabfromfile=F)
head(slist[[1]])

tr1 <- tabulateQ(qlist=slist)
tr1

#tabulateQ(tlist)
#tabulateQ(alist, writetable=TRUE)


sr1 <- summariseQ(tr1)
head(summariseQ(tabulateQ(slist)))

sfiles <- list.files(path=system.file("output_file_entire_ld.pruned_RUN1", package="pophelper"), full.names=T)


## load ind file 
inds <- read.delim(file ="structureindlabels.txt",header=FALSE,stringsAsFactors=F)
rownames(slist[[1]]) <- inds$V1

if(length(unique(sapply(slist,nrow)))==1) slist <- lapply(slist,"rownames<-",inds$V1)
lapply(slist, rownames)

## plotQ():create single-line barplots from qlist
p1 <- plotQ(slist[1],returnplot=T,exportplot=F,quiet=T,basesize=11,
            showindlab=T)  ###slist1[1]=> out_puts_fs_run1.1.png exported


# plot multiple runs separately
#p2<-plotQ(qlist=slist[1:3])


# join files into one figure
plotQ(qlist=readQ(sfiles)[c(1,2,3,4)], imgoutput="join")

# creating a short dataset
#slist1 <- sapply(slist,function(x) x[c(1:5,20:25,50:55,100:105,130:135),])
# normal usage
# p <- plotQ(slist1[1])
# modified for this document
#p <- plotQ(slist1[1],returnplot=T,exportplot=F,quiet=T,basesize=11)
#print(p$plot[[1]])

#Strip panel


#Cluster colours
# change colour of clusters
plotQ(qlist=slist[3:4], imgoutput="join", 
      clustercol=c("coral","steelblue","lightblue","purple","orange"))

#plotQ(qlist=slist1[16:17], imgoutput="join", 
#      clustercol=c("#A6CEE3", "#3F8EAA", "#79C360", "#E52829", "#FDB762","#ED8F47","#9471B4"))


clist <- list(
  "standard_12"=c("#2121D9","#9999FF","#DF0101","#04B404","#FFFB23","#FF9326","#A945FF","#0089B2","#B26314","#610B5E","#FE2E9A","#BFF217"),
  "rich.colors"=pophelper:::getColours(13))

# add length of palettes
lengths <- sapply(clist,length)

par(mar=c(0.2,4.5,0.2,0))
par(mfrow=c(length(clist),1))

for(i in 1:length(clist))
{
  {barplot(rep(1,max(lengths)),col=c(clist[[i]],rep("white",max(lengths)-length(clist[[i]]))),axes=F,border=F)
    text(x=-0.1,y=0.5,adj=1,label=names(clist)[i],xpd=T)}
}


clist <- list(
  "shiny"=c("#1D72F5","#DF0101","#77CE61", "#FF9326","#A945FF","#0089B2","#FDF060","#FFA6B2","#BFF217","#60D5FD","#CC1577","#F2B950","#7FB21D","#EC496F","#326397","#B26314","#027368","#A4A4A4","#610B5E"),
  "strong"=c("#11A4C8","#63C2C5","#1D4F9F","#0C516D","#2A2771","#396D35","#80C342","#725DA8","#B62025","#ED2224","#ED1943","#ED3995","#7E277C","#F7EC16","#F8941E","#8C2A1C","#808080"),
  "oceanfive"=c("#00A0B0", "#6A4A3C", "#CC333F", "#EB6841", "#EDC951"),
  "keeled"=c("#48B098", "#91CB62", "#FFEE3B", "#FB9013", "#FF3C28"),
  "vintage"=c("#400F13", "#027368", "#A3BF3F", "#F2B950", "#D93A2B"),
  "muted"=c("#46BDDD","#82DDCE","#F5F06A","#F5CC6A","#F57E6A"),
  "teal"=c("#CFF09E","#A8DBA8","#79BD9A","#3B8686","#0B486B"),
  "merry"=c("#5BC0EB","#FDE74C","#9BC53D","#E55934","#FA7921"),
  "funky"=c("#A6CEE3", "#3F8EAA", "#79C360", "#E52829", "#FDB762","#ED8F47","#9471B4"),
  "retro"=c("#01948E","#A9C4E2","#E23560","#01A7B3","#FDA963","#323665","#EC687D"),
  "cb_paired"=c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928"),
  "cb_set3"=c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#D9D9D9","#BC80BD","#CCEBC5","#FFED6F"),
  "morris"=c("#4D94CC","#34648A","#8B658A","#9ACD32","#CC95CC","#9ACD32","#8B3A39","#CD6601","#CC5C5B","#8A4500"),
  "wong"=c("#000000","#E69F00","#56B4E9","#009E73","#F0E442","#006699","#D55E00","#CC79A7"),
  "krzywinski"=c("#006E82","#8214A0","#005AC8","#00A0FA","#FA78FA","#14D2DC","#AA0A3C","#FA7850","#0AB45A","#F0F032","#A0FA82","#FAE6BE"))

# add length of palettes
lengths <- sapply(clist,length)
names(clist) <- paste0(names(clist),"_",lengths)


par(mar=c(0.2,6,0.2,0))
par(mfrow=c(length(clist),1))

for(i in 1:length(clist))
{
  {barplot(rep(1,max(lengths)),col=c(clist[[i]],rep("white",max(lengths)-length(clist[[i]]))),axes=F,border=F)
    text(x=-0.1,y=0.5,adj=1,label=names(clist)[i],xpd=T,cex=1.2)}
}


kelly <- c("#F2F3F4","#222222","#F3C300","#875692","#F38400","#A1CAF1","#BE0032","#C2B280","#848482","#008856","#E68FAC","#0067A5","#F99379","#604E97","#F6A600","#B3446C","#DCD300","#882D17","#8DB600", "#654522","#E25822","#2B3D26")

par(mar=c(0.2,4.5,0.2,0))
par(mfrow=c(2,1))

{barplot(rep(1,11),col=kelly[1:11],axes=F,border=F)
  text(x=-0.1,y=0.5,adj=1,label="kelly_22",xpd=T)
  barplot(rep(1,11),col=kelly[12:22],axes=F,border=F)}


p1 <- plotQ(slist[c(1,4:8)],imgoutput="join",returnplot=T,exportplot=F,quiet=T,basesize=11,
            clustercol=clist$shiny,splab=paste0("K=",sapply(slist[c(1,4:8)],ncol)))

p1$plot[[1]]

#Legend

# show legend
plotQ(qlist=slist[3], showlegend=T)
# move to right side
plotQ(qlist=slist1[3], showlegend=T, legendpos="right")
# change legend key size
plotQ(qlist=slist1[3], showlegend=T, legendkeysize=5)
# change legend text size
plotQ(qlist=slist1[3], showlegend=T, legendtextsize=5)



#Group labels

#plotQMultiline
p <- plotQMultiline(slist[4],exportplot=F,returnplot=T)
grid.arrange(p$plot[[1]][[1]])



