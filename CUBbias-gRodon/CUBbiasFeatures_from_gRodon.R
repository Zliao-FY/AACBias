# This module contains code(using gRodon) measuring CUB independent of length and composition(MILC) and CPB features

library(Biostrings)
library(gRodon)
library(dplyr)

path_to_genome_dir = 'path_to_your_bacteria_genome(.fna)_dir'
species_name <- list.files(path_to_genome_dir,pattern="*.fna")

for (file_name in species_name){
	path_to_speciesGenome_file = paste(path_to_genome_dir,file_name,sep='\\')
	genes <- readDNAStringSet(path_to_speciesGenome_file)

	rib_subunit_sum <- sum(grepl('ribosomal subunit',names(genes),ignore.case=T))
	rib_protein_sum <- sum(grepl('ribosomal protein',names(genes),ignore.case=T))
	
	if(rib_subunit_sum > rib_protein_sum){
		highly_expressed <- grepl('ribosomal subunit',names(genes),ignore.case=T)
	 } else{
		highly_expressed <- grepl('ribosomal protein',names(genes),ignore.case=T)
	 }
	
	species_predict <- predictGrowth(genes,highly_expressed)
	predict_grow <- data.frame(file_name,species_predict$CUBHE,species_predict$ConsistencyHE,species_predict$CPB,stringsAsFactors=FALSE)
	write.table(predict_grow, file = "organism_CUBbiasFeatures.txt", append = TRUE, sep = "\t",quote=FALSE, row.names = FALSE)
}