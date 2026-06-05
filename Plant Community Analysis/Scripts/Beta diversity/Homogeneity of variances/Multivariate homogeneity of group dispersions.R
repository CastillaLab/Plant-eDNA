################################################################################
# Script: Analysis of multivariate dispersion among lake plant communities
#
# Description:
# This script evaluates beta diversity dispersion (multivariate dispersion)
# among lake plant communities derived from eDNA metabarcoding data.
#
# The workflow:
#   1. Imports the plant community matrix
#   2. Filters out low-depth samples (< 1000 sequencing reads)
#   3. Calculates pairwise community dissimilarities among samples
#   4. Tests for homogeneity of multivariate dispersion among lakes
#      using PERMDISP (betadisper)
#   5. Performs ANOVA and permutation tests on dispersion
#   6. Conducts Tukey post hoc comparisons among lakes
#   7. Generates publication-quality visualizations of dispersion
#
# Required input file:
#   - FINAL_plants.median.clean.CM_0diff_NO_gen.rev_blast_28Dec21.csv
#
################################################################################


#######################################
# Load required packages
#######################################
library(ggplot2)   # Visualization
library(vegan)     # Community ecology and multivariate analysis
library(dplyr)     # Data manipulation


################################################################################
# 1. DATA PREPARATION
################################################################################

#######################################
# Import plant community matrix
#######################################

# Read metabarcoding community matrix
# Columns 1–6 contain metadata
# Columns 7–396 contain taxa read counts
CM<-read.csv("FINAL_plants.median.clean.CM_0diff_NO_gen.rev_blast_28Dec21.csv",header=TRUE)


################################################################################
# FILTER OUT SAMPLES WITH FEWER THAN 1000 READS
################################################################################

#######################################
# Calculate total reads per sample
#######################################

# Sum read counts across all taxa columns
CM$total=rowSums(CM[,7:396])

#######################################
# Explore sequencing depth distribution
#######################################

# Visualize read-depth distribution
hist(CM$total)

# Examine range of sequencing depths
range(CM$total)

#######################################
# Remove low-depth samples
#######################################

# Retain only samples with at least 1000 reads
CM.filtered<-filter(CM,total >= 1000)

#######################################
# Re-examine sequencing depth after filtering
#######################################

# Histogram of filtered samples
hist(CM.filtered$total)

# Range after filtering
range(CM.filtered$total)

################################################################################
# 2. ANALYSIS OF MULTIVARIATE DISPERSION (BETADISP)
################################################################################

#######################################
# Prepare community matrix
#######################################

# Retain:
#   - Column 1: lake identity
#   - Columns 7–396: taxa abundance data
CM.filtered<-CM.filtered[,c(1,7:396)]

# Check dimensions of filtered dataset
dim(CM.filtered)

#######################################
# Separate community matrix and groups
#######################################

# Community abundance matrix
# Each row represents a sample
# Each column represents a taxon
betadisp.community.all.taxa<-CM.filtered[,2:391]

# Grouping variable defining lake identity
betadips.groups.all.taxa<-CM.filtered$lake_name

#######################################
# Calculate ecological dissimilarities
#######################################

# Compute pairwise Bray–Curtis dissimilarities
# among all samples using vegan::vegdist()
#
# Bray–Curtis distances quantify compositional
# differences among plant communities based on
# taxa abundances.
dis.betadisp.all.taxa<-vegdist(betadisp.community.all.taxa)


#######################################
# Evaluate homogeneity of dispersion
#######################################

# vegan::betadisper() calculates the distance
# of each sample to the centroid (spatial median)
# of its corresponding lake.
#
# Larger distances indicate greater within-lake
# heterogeneity in plant community composition.
#
# sqrt.dist = TRUE applies square-root
# transformation to distances to improve
# ordination properties.
mod.betadisp.all.taxa<-betadisper(dis.betadisp.all.taxa,betadips.groups.all.taxa,sqrt.dist = TRUE)

# Display model summary
mod.betadisp.all.taxa

#######################################
# Test significance of dispersion differences
#######################################

# Parametric ANOVA test
anova(mod.betadisp.all.taxa)

# Permutation-based test using 999 permutations
#
# pairwise = TRUE performs pairwise comparisons
# among lakes.
permutest(mod.betadisp.all.taxa, pairwise=TRUE, permutations=999)


#######################################
# Post hoc pairwise comparisons
#######################################

# Tukey Honest Significant Difference test
# compares dispersion among all lake pairs
(mod.HSD<-TukeyHSD(mod.betadisp.all.taxa))

# Plot confidence intervals from Tukey test
plot(mod.HSD, yaxt="n")

#######################################
# Ordination plot of multivariate dispersion
#######################################

# Visualize samples and lake centroids
# on the first two principal coordinates axes
plot(mod.betadisp.all.taxa, main="Dispersion across lakes", label=FALSE)

################################################################################
# 3. VISUALIZATION OF DISTANCES TO GROUP CENTROIDS
################################################################################

#######################################
# Create dataframe of distances
#######################################

# Extract:
#   - Lake identity
#   - Distance of each sample to its lake centroid
df.dist.cent.all.taxa<-data.frame(mod.betadisp.all.taxa$group,
                                  mod.betadisp.all.taxa$distances)


#######################################
# Generate boxplot of dispersion
#######################################
p1<-ggplot(df.dist.cent.all.taxa,aes(x = gsub("_", " ", mod.betadisp.all.taxa.group),
                                 y=mod.betadisp.all.taxa.distances))+
  stat_boxplot(geom = "errorbar")+geom_boxplot()+theme_classic()+
  ylab("Distance to spatial median")+
  ylim(0,1)+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 60,hjust = 1, size=12))


#################
# Save figure
#################
ggsave("betadisp_distances_plot.png", p1, width = 6, height = 5, dpi = 600)
