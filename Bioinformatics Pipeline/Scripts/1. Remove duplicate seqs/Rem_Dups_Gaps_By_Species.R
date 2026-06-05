#06/21/19  Script to take a fasta file, calculate nucleotide distances within species and keep only the longest 
#of duplicate seqs in the output alignment (or one of the dups if they are all the same length)

#Load packages
library(ape)
library(haplotypes)

setwd("C:/JK/2018 eDNA metabarcoding marker dev- John R/rbcL - Plants/Sequences/Remove Dups folder")

#assign the infile
infile = "Plant_rbcL_align_PrimerTrunc_23_081919.fas"

#create the outfile
outfile = "rbcL_nodup_outfile_081919.fas"

#read in the infile (changes N's to ?'s and gives you a warning)
dat <- read.fas(infile)

#isolate the sequence ID (Acc#_Genus_species)
sppnames <- names(dat)

#separate the ID to get genus and species and combine into one string separated by _
getname <- function(x) {
	split.name <- strsplit(x, split = "_")[[1]]
	paste.name <- paste0(split.name[2],"_",split.name[3])
	paste.name
}

spp <- sapply(sppnames, FUN= getname, USE.NAMES = FALSE)

out.df <- data.frame(
	seqID = c(1:length(spp)),
	species = spp, 
	seqlen = rowSums(dat@sequence != "?" & dat@sequence != "-"), 
	keep = rep(TRUE, length(spp))
)

#Loop over the unique species in the alignment
for(sp in unique(spp)) {
	
	#Extract just the portion of out.df that specifically corresponds to this species
	tmp.df <- out.df[out.df$species == sp,]
	
	#Only look for duplicates if there's more than 1 sequence in the alignment
	if(length(tmp.df$species) > 1) {

	#Distance matrix calculation using gaps as 5th character, make into matrix format
	#and keep only the bottom half of the matrix
		dmat <- distance(as.dna(dat[spp == sp,]), indels = "5")
		dmat1 <- as.matrix(dmat)
		dmat1[upper.tri(dmat1, diag = TRUE)] <- NA
		
	#Are there duplicates of a given sequence?
		is.dup <- colSums(dmat1 == 0, na.rm = TRUE)
		
	#If so, which sequences are duplicated?  Which is longest?  Keep only the longest!
		for(x in 1:length(is.dup)) {
			if(is.dup[x] > 0) {
				dupl <- c(x,which(dmat1[,x] == 0))
				duplID <- tmp.df$seqID[dupl]
				tmp2.df <- tmp.df[tmp.df$seqID %in% duplID,]
				keeper <- tmp2.df$seqID[tmp2.df$seqlen == max(tmp2.df$seqlen)]
				duplID <- duplID[duplID != keeper[1]]
				out.df$keep[out.df$seqID %in% duplID] = FALSE
				rm(tmp2.df, dupl, duplID, keeper)
			}
		}

		rm(dmat, dmat1, is.dup)	
	
	} else {
		rm(tmp.df)
	}
}

#make an object containing all those seqs to keep
nodup_dat <- dat[which(out.df$keep == TRUE),]

#write the screened alignment to file
write.dna(nodup_dat, file = outfile, format = "fasta")

#the end!