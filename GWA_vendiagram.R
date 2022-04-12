

require(VennDiagram)

##### Worked
### source: https://stackoverflow.com/questions/24736637/r-color-overlaps-in-venn-diagram-by-size-of-overlap/26821445



####
setwd("~/Documents/malaria_vector_reserach/input_data/a.aegypti/ven_diagram_all_scans")
vars<-c("latitude_e","longitude_e","MAT_e","MWMT_e","MCMT_e","TD_e","MAP_e","MSP_e","AHM_e","SHM_e","DD_0_e","DD5_e","DD_18_e","DD18_e","NFFD_e","bFFP_e","eFFP_e","FFP_e","PAS_e","EMT_e","EXT_e","Eref_e","CMD_e","MAR_e","RH_e")

for(i in 1:length(vars)){
  
    
  baypass<-read.table(file = paste("input_ven_top_candidate_loci_BayPass/baypass_A.aegypi_12pops_BF20_list_candidate_loci_",vars[i],".txt", sep = ""),header = TRUE)

    lfmm<-read.table(file = paste("input_ven_top_candidate_loci_lfmm_FDR_0.01/lfmm_A.aegypti_12_pops_list_candidate_loci_FDR_0.01_",vars[i],".txt", sep = ""),header = TRUE)

    xtx<-read.table(file = "baypass_XtX_A.aegypi_12pops_POD22_list_candidate_loci.txt", header = FALSE)
    colnames(xtx)<-c("snp_id","CHR","pos","xtx")

    pcadapt<-read.table(file = "pcadapt_A.aegypti_12_pops_list_candidate_loci_FDR_0.01_from_whole_genome.txt", header = TRUE)


   baypass_ven<-baypass$snp_id
   lfmm_ven<-lfmm$snp_id
   xtx_ven<-xtx$snp_id
   pcadapt_ven<-pcadapt$snp_id


   pdf(file = paste("vediagram_all_GWS_A.aegypti_12_pops_",vars[i],".pdf",sep=""))
   v1 <- venn.diagram(list(BayPass=baypass_ven, LFMM=lfmm_ven, XtX=xtx_ven, PcAdapt=pcadapt_ven), filename=NULL, fill=rainbow(4))
   grid.newpage()
   grid.draw(v1)
   dev.off()

}




