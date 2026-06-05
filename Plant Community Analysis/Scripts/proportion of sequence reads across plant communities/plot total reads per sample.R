################################################################################
# Script: Read depth visualization for eDNA plant metabarcoding data
# 
# This script imports a metabarcoding data set and species metadata, filters 
# the data set to retain only taxa identified at species level, calculates the
# total number of sequencing reads per sample, and generate two plots:
#
# 1. Histogram of total reads per sample
# 2. Boxplot of total reads per lake
#
# Required input files:
#
# - FINAL_plants.median.clean.CM_0diff_NO_gen.rev_blast_28Dec21.csv: 
#   Plant community matrix with the reads for each taxa detected in the eDNA
#   metabarcoding analysis
#
# - spinfo.csv: Metadata with the list of species detected to only use
#   detections at the species level.
#
################################################################################

#######################################
# Load required packages
#######################################

library(ggplot2)    # Data visualization
library(cowplot)    # Combining plots
library(dplyr)      # Data manipulation


#######################################
# Import input data
#######################################

# Read metabarcoding dataset
df1<-read.csv("FINAL_plants.median.clean.CM_0diff_NO_gen.rev_blast_28Dec21.csv", header=TRUE)

# Read species information table
spinfo <- read.csv("spinfo.csv")


#########################################################################
# Filter the dataset to retain only taxa corresponding to species-level 
# identifications.
#########################################################################

# Extract species names from metadata table
keep_spp<-spinfo$species

# Keep:
#   - First six metadata columns
#   - Species columns matching the target species list
df1<-df1[,c(1:6,which(colnames(df1) %in% keep_spp))]


#######################################
# Calculate total reads per sample
#######################################

# Sum read counts across all species columns
# Species data start at column 7
df1$total.reads = rowSums(df1[,7:227])


#######################################
# Format variables
#######################################

# Replace underscores with spaces in lake names
df1$lake_name <- gsub("_", " ", df1$lake_name)

# Convert categorical variables to factors
df1$lake_name<-as.factor(df1$lake_name)
df1$sample<-as.factor(df1$sample)
df1$sample_type<-as.factor(df1$sample_type)

#######################################
# Plot 1: Histogram of total reads
#######################################
p1<-ggplot(df1,aes(x=total.reads))+
  geom_histogram(binwidth = 500, color="black",
                 fill="white")+
  theme_classic() + xlab("Total number of reads")+
  ylim(0,250)+
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size = 16))+
  geom_vline(xintercept = 1000, linetype="dashed",
             color="red", size=1)
  

#######################################
# Plot 2: Total reads by lake
#######################################
p2<-ggplot(df1,aes(x=lake_name,y=total.reads))+
  geom_boxplot()+theme_classic()+ theme(axis.text.x = element_text(angle = 60,hjust = 1, size=8),
                                        axis.text.y = element_text(size=14),
                                        axis.title = element_text(size=16))+
  ylab("Total number \n of sequence reads")+
  xlab("Lake")+
  ylim(0,30000)


#######################################
# Combine plots into a single figure
#######################################
combinedPlot<-plot_grid(p1,p2,ncol=2,nrow=1,labels='auto',label_size = 20,
                        align='hv',rel_widths=1,rel_heights=1)



#######################################
# Export figure
#######################################
png("number of reads per sample.png",width = 4000, height = 1600, units = "px", res=300)
combinedPlot
dev.off()


