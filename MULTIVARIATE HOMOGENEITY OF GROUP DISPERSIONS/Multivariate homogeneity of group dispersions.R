# First, load all packages
library(ggplot2)
library(cowplot)
library(dplyr)
library(vegan)

############################
##   1. DATA PREPARATION  ##  
############################

# Load the plant community matrix
CM<-read.csv("FINAL_plants.median.clean.CM_0diff_NO_gen.rev_blast_28Dec21.csv",header=TRUE)

######################################################
#  FILTERING OUT SAMPLES WITH FEWER THAN 1000 READS  #
######################################################

# Create a new variable total with the total number of reads per sample
CM$total=rowSums(CM[,4:393])

# Let's check the frequency distribution and range of total
hist(CM$total)
range(CM$total)

# Filter the plant comunity matrix to keep only samples with >= 1000 sequences
CM.filtered<-filter(CM,total >= 1000)

# Let's check the frequency distribution and range of total
hist(CM.filtered$total)
range(CM.filtered$total)

#############################################################
# Multivariate homogeneity of group dispersions (variances) #
#############################################################

#############################################
#  All taxa - including above-species taxa  #
#############################################

# Keep only lake ID and the counts for each taxa
CM.filt.counts<-CM.filtered[,c(1,4:393)]
dim(CM.filt.counts)

# Remove the lake ID
betadisp.community.all.taxa<-CM.filt.counts[,2:391]
betadips.groups.all.taxa<-CM.filt.counts$lake_name

dis.betadisp.all.taxa<-vegdist(betadisp.community.all.taxa)

#Examine whether lakes exhibit homogeneity in their dispersions (variances)
mod.betadisp.all.taxa<-betadisper(dis.betadisp.all.taxa,betadips.groups.all.taxa,sqrt.dist = TRUE)
mod.betadisp.all.taxa

#Perform anova test
anova(mod.betadisp.all.taxa)

#Permutation test for F
permutest(mod.betadisp.all.taxa, pairwise=TRUE, permutations=999)

#Tukey's Honest Significance Test
(mod.HSD<-TukeyHSD(mod.betadisp.all.taxa))
plot(mod.HSD, yaxt="n")

#Plot the groups and distances to the centroids on the first two PCoA axes
plot(mod.betadisp.all.taxa, main="Dispersion across lakes", label=FALSE)

#Create a dataframe with the group and the distances to the spatial median

df.dist.cent.all.taxa<-data.frame(mod.betadisp.all.taxa$group,mod.betadisp.all.taxa$distances)

p1<-ggplot(df.dist.cent.all.taxa,aes(x=mod.betadisp.all.taxa.group,
                                                         y=mod.betadisp.all.taxa.distances))+
  stat_boxplot(geom = "errorbar")+geom_boxplot()+theme_classic()+
  ylab("Distance to spatial median")+
  ylim(0,1)+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 60,hjust = 1, size=12))

png("Non-euclidean distances plot.png",width = 4000, height = 3200, units = "px", res=300)
p1
dev.off()
