
#Install DECIPHER if necessary

#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
#BiocManager::install("DECIPHER")

#Load DECIPHER

library(DECIPHER)


#Help while in DECIPHER is:  ? Function/term to look up


#Set working directory (put your alignment in here)

setwd("C:/JK/2018 eDNA metabarcoding marker dev- John R/rbcL - Plants/Sequences/R_align")


#set the Genetic Code (see https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi for a description of the different codes, and their #'s)

x <- getGeneticCode(id_or_name2="11", full.search=FALSE, as.data.frame=FALSE)


#Align protein coding sequences using a single genetic code (Your alignment must have all gap (and .) characters removed both before and after DNA sequences)

#Read in the alignment

dna <- readDNAStringSet("align.plant.rbcl.seqs_042519_LP_JKnogaps1.fas")

#Call the AlignTranslation function with your desired arguments

DNA <- AlignTranslation(dna, sense = "+", direction = "5' to 3'", readingFrame = NA, type = "DNAStringSet", geneticCode = x)

#write the outfile to your working directory

writeXStringSet(DNA, file = "DNA_out.fas")

#Done
