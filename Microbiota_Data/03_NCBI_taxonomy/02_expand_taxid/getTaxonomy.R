#https://cran.r-project.org/web/packages/taxonomizr/readme/README.html

#install.packages("taxonomizr")
library(taxonomizr) #packageVersion("taxonomizr")

## set source to file location
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))

#clear the object from memory
rm(list=ls())

#install DB
#taxonomizr::prepareDatabase('accessionTaxa.sql')


##get taxonomy based on taxID (taxID was blasted before)
#we allowed max 5 hits per query
#during blast we only allowed very accurated hits with evalue 1e-5 or higher
taxIDs <- read.delim("../01_blast_taxid/taxid.tab", sep = "\t", header = F)
TAXA_NCBI_all <- taxonomizr::getTaxonomy(taxIDs[,3],'taxonomizr/accessionTaxa.sql')
TAXA_NCBI_all <- as.data.frame(TAXA_NCBI_all)
TAXA_NCBI_all$OTU <- taxIDs$V1

#selected for the hit with the most entries
OTUs <- unique(taxIDs[,1])
levels <- colnames(TAXA_NCBI_all)
TAXA <- data.frame(matrix(nrow=length(OTUs), ncol=length(levels),dimnames = list(OTUs, levels)))
for(OTU in OTUs){
  TAXA_NCBI_sub <- TAXA_NCBI_all[TAXA_NCBI_all$OTU == OTU,]
  #get the taxa with most entry (if  multiple with equal entry, take the first one)
  nNA_perOTU <- rowSums(is.na(TAXA_NCBI_sub))
  hit_selected <- which(nNA_perOTU == min(nNA_perOTU))[1]
  #taxa of seleceted hit
  TAXA[rownames(TAXA)==OTU,] <- TAXA_NCBI_sub[hit_selected,]
}

#add classifier and remove OTU col
TAXA <- TAXA[,colnames(TAXA)!="OTU"]
TAXA$classifier <- "NCBI"

#rename superkingdom to kingdom
colnames(TAXA) <- gsub("superkingdom", "kingdom", colnames(TAXA))

write.table(TAXA, "../bacteria_TAXA100_NCBI.tab", sep="\t")

