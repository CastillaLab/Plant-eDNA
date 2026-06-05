#originally created on: August 27, 2019
#by: Jeannette Kanefsky

#ABOUT: This script was written to extract taxonomies out of the ITIS database (www.ITIS.gov) for a list of 
#species in an alignment and create a taxonomy file for Plants for Mothur.  It is designed to
#insert NA where taxonomies can not be retrieved; they will have to be edited manually in the
#final taxonomy file.

#It can be edited to adapt to different animal groups by altering accepted size of classification output and
#changing what taxonomic levels are saved (animals have phyla, plants have divisions).  ***You need to go 
#to ITIS and see what taxonomic levels are reported for your group of interest.***

#loading in libraries

library(dplyr)
library(tidyr)
library(stringr)
library(taxize)


#setting working directory
setwd("C:/Users/kanefsk1/Desktop/R Taxize taxonomy prep scripting")

#Listing files (this script requires Input and Output folders within the working directory)
list.files("Input/")

#Read in the sequence alignment and create a dataframe
df <- readLines("Input/test_rbcL_082219.fas")
df <- df[grepl(pattern = ">", x = df)]
df <- data.frame(ids = df, stringsAsFactors = F)

df$ids2 <- df$ids
head(df)

#Split the sequence names into Accesion #s, Genus and Species names and put Genus and Species together to make scientific names 
df <- separate(data = df, col = ids2, into = c("acc","genus","species"), sep = "_")
df$sci.name <- paste(df$genus, df$species)
df$accNo <- gsub(pattern = ">", replacement = "", x = df$acc)
head(df)


#initialize an empty dataframe to collect all of the seq names and their taxonomies
TaxFile <- data.frame(AccNO=character(), Taxonomy=character(), stringsAsFactors = FALSE)

#create an empty Taxa variable
Taxa <- ""

#create empty variables to hold the taxonomic names (tailor to your group)
kingdom <- ""
division <- ""
class <- ""
order <- ""
family <- ""
genus <- ""
species <- ""

#make a fake placeholder dataframe for tax.temp1 if the species is not found in ITIS and the classification file is NULL
fake.df <- data.frame(AccNO=character(), Taxonomy=character(), stringsAsFactors = FALSE)
fake.df[nrow(fake.df) + 1,] <- c("Empty", "ReallyEmpty")

#Iterate through the species list, save seq names and species names to variables and send to the ITIS website database and retrieve the classification
### NOTE:  You will have to watch the output during the retrieval of the classifications from ITIS because if more than one match
#is made for a name (for example, when there are varieties or subspecies, or for some reason the same name with 2 different TSN #'s), the 
#process will pause for you to input the row # of the taxon you want a classification retrieved for (usually the one whose name usage is accepted)

#for plants (edit for your group)
for (row in 1:nrow(df)) {
  seq.name <- as.character(gsub(pattern = ">", replacement = "", x = df[row, "ids"]))  
  species.name <- as.character(df[row, "sci.name"])
  #print(seq.name)
  #print(species.name)
  tax.temp <- classification(species.name, db = 'itis')
  tax.temp1 <- rbind(tax.temp)
  if (is.null(tax.temp1)) {
    tax.temp1 <- fake.df
    }
  #print(tax.temp1)
  
  #for valid classification files, assign the names at the appropriate ranks to variables  
  if (nrow(tax.temp1) > 1) {
    kingdom <- as.character(tax.temp1 %>% filter(rank == "kingdom"))[1]
    division <- as.character(tax.temp1 %>% filter(rank == "division"))[1]
    class <- as.character(tax.temp1 %>% filter(rank == "class"))[1]
    order <- as.character(tax.temp1 %>% filter(rank == "order"))[1]
    family <- as.character(tax.temp1 %>% filter(rank == "family"))[1]
    genus <- as.character(tax.temp1 %>% filter(rank == "genus"))[1]
    species <- as.character(tax.temp1 %>% filter(rank == "species"))[1]
    species.only <- unlist(strsplit(species, " "))[2]
    
  }
  
  #for classification files for species that don't have valid names, use NA for all levels to denote they have to researched by hand  
  else {
    kingdom <- "NA"
    division <- "NA"
    class <- "NA"
    order <- "NA"
    family <- "NA"
    genus <- "NA"
    species.only <- "NA"
  }
  
  #concatenate the names at all ranks with ;'s and assign to a variable, then add a row containing the seq_species name and 
  #classification for that species to the TaxFile
  Taxa <- as.character(paste0(kingdom, ";", division, ";", class, ";", order, ";", family, ";", genus, ";", species.only, ";"))
  #print(Taxa)
  TaxFile[nrow(TaxFile) + 1,] <- c(seq.name, Taxa)
   
  #reset taxonomic variables to empty (edit for your group)
  kingdom <- ""
  division <- ""
  class <- ""
  order <- ""
  family <- ""
  genus <- ""
  species <- ""
  
  #reset Taxa variable to empty
  Taxa <- ""
  
  #reset tax.temp, tax.temp1, seq.name and species.name to empty
  tax.temp <- ""
  tax.temp1 <- ""
  seq.name <- ""
  species.name <- ""
  
}

#writing the full taxonomy (TaxFile) to tab delimited text file (edit filename for group and gene name)
write.table(x = TaxFile, file = "Output/PlantrbcL.new.taxonomy.file.txt", append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

#Done!