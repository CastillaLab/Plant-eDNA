# Run FindNCBI function to retrieve (plant rbcL) sequences from GenBank


#load in functions
source("Find_NCBI.R")

#read in species list for FindNCBI
sl <- read.table("Input/All_Plant_rbcL_list_051819.txt", header = T, stringsAsFactors = F)
sl <- data.frame(sci.name = paste(sl$sample_number,sl$sci.name,sep = " "))

#run FindNCBI (set locus name and maximum # of seqs to retrieve per species)
FindNCBI(species_list = sl, locus = "rcbL", setmax = 10)
