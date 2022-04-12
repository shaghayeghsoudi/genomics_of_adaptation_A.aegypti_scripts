setwd("~/Documents/malaria_vector_reserach/input_data/a.aegypti/Go-enrichment/intersect_genes_three_ways_genomescans")
library("biomaRt")
library(topGO)

#https://datacatz.wordpress.com/2018/01/19/gene-set-enrichment-analysis-with-topgo-part-1/
#https://www.biostars.org/p/1226/ 
#https://www.biostars.org/p/250927/
#https://www.biostars.org/p/247303/

  
  #collect gene names from biomart
#mart <- biomaRt::useMart(biomart = "plants_mart",
#                         dataset = "athaliana_eg_gene",
#                         host = 'plants.ensembl.org')


#listDatasets(wormbase)

mart <- useMart(biomart = "metazoa_mart", 
                    dataset = "aalvpagwg_eg_gene",
                    host = "https://metazoa.ensembl.org", 
                    port = 443)

mart <- biomaRt::useMart(biomart = "metazoa_mart",
                         dataset = "aalvpagwg_eg_gene",
                         host = 'https://metazoa.ensembl.org')

#listAttributes(mart) ## ensembl_gene_id and go_id selected from getBM
GTOGO <- getBM(attributes = c( "ensembl_gene_id",
                                        "go_id"), mart = mart)


#examine result
head (GTOGO)
#Remove blank entries
GTOGO <- GTOGO[GTOGO$go_id != '',]
# convert from table format to list format
geneID2GO <- by(GTOGO$go_id,
                GTOGO$ensembl_gene_id,
                function(x) as.character(x))
#examine result
head (geneID2GO)


#Make topGO data object
all.genes <- sort(unique(as.character(GTOGO$ensembl_gene_id)))
#int.genes <- sample(x = all.genes, size = 200) # some random genes 

intersect_list<-read.table(file = "overlapped_genesALL_VARS_UNIQUE.txt", header  = FALSE)
int.genes<-intersect_list[,5]

int.genes <- factor(as.integer(all.genes %in% int.genes))
names(int.genes) = all.genes

go.obj <- new("topGOdata", ontology='BP'
              , allGenes = int.genes
              , annot = annFUN.gene2GO
              , gene2GO = geneID2GO
)

results <- runTest(go.obj, algorithm = "elim", statistic = "fisher")
results.tab <- GenTable(object = go.obj, elimFisher = results)
# Plot results using the GO hierarchical DAG:
showSigOfNodes(go.obj, score(results), firstSigNode=15, useInfo ='all')


#####

GTOGO
bm <- GTOGO[!duplicated(GTOGO[,1]),]
nrow(bm)

#Extract the gene-level statistics that we will use:
d_piano <- merge(diffExpRes, bm, by="ensembl_gene_id", all.x=T, sort=F)



