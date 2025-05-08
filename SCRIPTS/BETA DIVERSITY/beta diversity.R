# First, load all packages
library(ggplot2)
library(cowplot)
library(dplyr)
library(vegan)
library(betapart)


#library(tidyverse)

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

####################################################
# SUBSET TERRESTRIAL, AQUATIC, AND WETLAND PLANTS  #
####################################################

# Load the list with the species, their life habit, and status
spinfo<-read.csv("spinfo.csv", header=TRUE)

# We create a list with all species 
keep_spp_entire<-spinfo$species[spinfo$community == "entire"]
length(keep_spp_entire)

# We keep only all species in our plant community matrix
CM.filt.entire<-CM.filtered[,c(1:3,which(colnames(CM.filtered) %in% keep_spp_entire))]
head(CM.filt.entire)

# We create a list with the terrestrial species 
keep_spp_terr<-spinfo$species[spinfo$habitat == "terrestrial"]
length(keep_spp_terr)

# We keep only the terrestrial species in our plant community matrix
CM.filt.terr<-CM.filtered[,c(1:3,which(colnames(CM.filtered) %in% keep_spp_terr))]
head(CM.filt.terr)


# We create a list with the wetland species 
keep_spp_wet<-spinfo$species[spinfo$habitat == "wetland"]
length(keep_spp_wet)

# We keep only the wetland species in our plant community matrix
CM.filt.wet<-CM.filtered[,c(1:3,which(colnames(CM.filtered) %in% keep_spp_wet))]
head(CM.filt.wet)

# We create a list with the aquatic species 
keep_spp_aqu<-spinfo$species[spinfo$habitat == "aquatic"]
length(keep_spp_aqu)

# We keep only the wetland species in our plant community matrix
CM.filt.aqu<-CM.filtered[,c(1:3,which(colnames(CM.filtered) %in% keep_spp_aqu))]
head(CM.filt.aqu)

######################################
##     Summarize for each lake       #
######################################

# Entire community
CM.filt.entire<-CM.filt.entire[,c(1,4:222)]
CM.filt.entire.lake<-CM.filt.entire %>% group_by(lake_name) %>% summarise(across(everything(),list(sum)))

write.csv(CM.filt.entire.lake,"CM.filt.entire.lake.csv", row.names=FALSE)


# Terrestrial plants
CM.filt.terr<-CM.filt.terr[,c(1,4:132)]
CM.filt.terr.lake<-CM.filt.terr %>% group_by(lake_name) %>% summarise(across(everything(),list(sum)))

write.csv(CM.filt.terr.lake,"CM.filt.terr.lake.csv", row.names=FALSE)

# Wetland plants
CM.filt.wet<-CM.filt.wet[,c(1,4:63)]
CM.filt.wet.lake<-CM.filt.wet %>% group_by(lake_name) %>% summarise(across(everything(),list(sum)))

write.csv(CM.filt.wet.lake,"CM.filt.wet.lake.csv", row.names=FALSE)

# Aquatic plants
CM.filt.aqu<-CM.filt.aqu[,c(1,4:33)]
CM.filt.aqu.lake<-CM.filt.aqu %>% group_by(lake_name) %>% summarise(across(everything(),list(sum)))

write.csv(CM.filt.aqu.lake,"CM.filt.aqu.lake.csv", row.names=FALSE)


###########################################
##      2. PERMANOVA AND ddRAD           ##  
###########################################

##################################
#      2.1 Entire community      #  
##################################

CM.filt.entire.lake<-read.csv("CM.filt.entire.lake.csv", header=TRUE)
df.entire<-CM.filt.entire.lake[,2:220]

lake.env<-read.csv("eDNA.env.csv",header=TRUE)

lake.entire.mdf <- as.matrix.data.frame(df.entire)
rownames(lake.entire.mdf) <- lake.env$lake.name

lake.bray.entire <- vegdist(df.entire, method = "bray")
lake.bray.entire

lake.jac.entire <- vegdist(df.entire, method = "jaccard", binary = T)
lake.jac.entire

## Principal Coordinates Analysis (PCoA)  

# calculate principal coordinates analysis (Bray-Curtis)
pcoa.entire.bray <- cmdscale(lake.bray.entire, k = 2, eig = T)

# extract axis positions for each site from cmdscale object and create a dataframe for plotting
pcoa.entire.bray.plotting <- as.data.frame(pcoa.entire.bray$points)
colnames(pcoa.entire.bray.plotting) <- c("axis_1", "axis_2")
pcoa.entire.bray.plotting$site <- rownames(pcoa.entire.bray.plotting)

# calculate the proportion of variance in the data which is explained by the first two PCoA axes
pcoa.entire.bray$eig[1]/(sum(pcoa.entire.bray$eig))
pcoa.entire.bray$eig[2]/(sum(pcoa.entire.bray$eig))

# create a PCoA plot
pcoa.entire.bray.plot <- ggplot(pcoa.entire.bray.plotting, aes(x = axis_1, y = axis_2, colour = site)) +
  geom_point(size = 3) +
  scale_colour_viridis_d(option = "magma", begin = 0.2, end = 0.8) +
  theme_bw() + 
  xlab("PCoA 1 (16.6%)") +
  ylab("PCoA 2 (15.9%)") +
  annotate(geom = 'text', label = 'Bray-Curtis', x = Inf, y = -Inf, hjust = 1.15, vjust = -1)

# repeat process with Jaccard dissimilarity matrix
pcoa.entire.jac <- cmdscale(lake.jac.entire, k = 2, eig = T)

pcoa.entire.jac.plotting <- as.data.frame(pcoa.entire.jac$points)
colnames(pcoa.entire.jac.plotting) <- c("axis_1", "axis_2")
pcoa.entire.jac.plotting$site <- rownames(pcoa.entire.jac.plotting)

pcoa.entire.jac$eig[1]/(sum(pcoa.entire.jac$eig))
pcoa.entire.jac$eig[2]/(sum(pcoa.entire.jac$eig))

pcoa.entire.jac.plot <- ggplot(pcoa.entire.jac.plotting, aes(x = axis_1, y = axis_2, colour = site)) +
  geom_point(size = 3) +
  scale_colour_viridis_d(option = "magma", begin = 0.2, end = 0.8) +
  theme_bw() + 
  xlab("PCoA 1 (15.5%)") +
  ylab("PCoA 2 (10.8%)") +
  annotate(geom = 'text', label = 'Jaccard', x = Inf, y = -Inf, hjust = 1.215, vjust = -1)

# extract plot legend
legend <- get_legend(pcoa.entire.jac.plot)

# plot Bray-Curtis PCoA and Jaccard PCoA side by side
pcoa.bray.jac.entire<-plot_grid(pcoa.entire.bray.plot + theme(legend.position = 'none'), pcoa.entire.jac.plot + theme(legend.position = 'none'), legend, ncol = 3, rel_widths = c(1,1,0.5))

pcoa.bray.jac.entire

##  Permanova

entire.permanova<-adonis2(lake.bray.entire ~ nhd_lat + julian.day + lake.area + IWS.stream.density + area.upstream.lakes + perc.agric.dev.land + max.depth,
       data = lake.env, permutations = 999)

##  ddRDA

lake.env.z <- lake.env

lake.env.z$julian.day <- (lake.env.z$julian.day - mean(lake.env.z$julian.day))/sd(lake.env.z$julian.day)
lake.env.z$max.depth <- (lake.env.z$max.depth - mean(lake.env.z$max.depth))/sd(lake.env.z$max.depth)
lake.env.z$lake.area <- (lake.env.z$lake.area  - mean(lake.env.z$lake.area))/sd(lake.env.z$lake.area)
lake.env.z$IWS.stream.density <- (lake.env.z$IWS.stream.density - mean(lake.env.z$IWS.stream.density))/sd(lake.env.z$IWS.stream.density)
lake.env.z$area.upstream.lakes <- (lake.env.z$area.upstream.lakes - mean(lake.env.z$area.upstream.lakes))/sd(lake.env.z$area.upstream.lakes)
lake.env.z$perc.agric.dev.land <- (lake.env.z$perc.agric.dev.land - mean(lake.env.z$perc.agric.dev.land))/sd(lake.env.z$perc.agric.dev.land)
lake.env.z$nhd_lat <- (lake.env.z$nhd_lat - mean(lake.env.z$nhd_lat))/sd(lake.env.z$nhd_lat)

dbRDA.entire.full <- capscale(lake.bray.entire ~ nhd_lat + julian.day+lake.area+IWS.stream.density+area.upstream.lakes+perc.agric.dev.land+max.depth,
                              lake.env.z)

vif.cca(dbRDA.entire.full)

anova(dbRDA.entire.full, by = "terms")
summary(dbRDA.entire.full)

smry.entire <- summary(dbRDA.entire.full)
scrs.entire <- scores(dbRDA.entire.full)
df.ddRAD.entire  <- data.frame(smry.entire$sites[,1:2]) # site scores for RDA1 and RDA2
df.ddRAD.entire$site <- rownames(df.ddRAD.entire)  #add site names
df.ddRAD.entire$site <-as.factor(df.ddRAD.entire$site)  #add site names

df.ddRAD.entire2  <- data.frame(smry.entire$biplot[,1:2])  # mapping environmental variables

rda.entire.plot <- ggplot(df.ddRAD.entire, aes(x=CAP1, y=CAP2, colour = site)) + 
  geom_segment(data=df.ddRAD.entire2, aes(x=0, xend=CAP1, y=0, yend=CAP2), 
               color="grey50", arrow=arrow(length=unit(0.01,"npc"))) +
  geom_text(data=df.ddRAD.entire2, aes(x=CAP1,y=CAP2,label=rownames(df.ddRAD.entire2),
                                       hjust=0.5*(1-sign(CAP1)),vjust=0.5*(1-sign(CAP2))), 
            color="grey50", size=3) +
  geom_point(size = 3) +
  scale_colour_viridis_d(option = "magma", begin = 0.2, end = 0.8) +
  geom_text(aes(label=rownames(df.ddRAD.entire),
                hjust=0,vjust=1.5), colour = "black",size=3) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  xlim(-2, 2) +
  ylim(-2, 1.5) +
  xlab("RDA1 (33.5%)") + # this percentage comes from the CAP1 'importance of components:' proportion explained, which can be found in summary(dbRDA.mat) 
  ylab("RDA2 (23.7%)") + # this percentage comes from the CAP2 'importance of components:' proportion explained, which can be found in summary(dbRDA.mat) 
  coord_fixed() +
  theme_bw() 
rda.entire.plot

#######################################
#      2.2 Terrestrial community      #  
#######################################

CM.filt.terr.lake<-read.csv("CM.filt.terr.lake.csv", header=TRUE)

df.terr<-CM.filt.terr.lake[,2:130]

lake.terr.mdf <- as.matrix.data.frame(df.terr)
rownames(lake.terr.mdf) <- lake.env$lake.name

lake.bray.terr <- vegdist(df.terr, method = "bray")
lake.bray.terr

lake.jac.terr <- vegdist(df.terr, method = "jaccard", binary = T)
lake.jac.terr

## Principal Coordinates Analysis (PCoA)

# calculate principal coordinates analysis (Bray-Curtis)
pcoa.terr.bray <- cmdscale(lake.bray.terr, k = 2, eig = T)

# extract axis positions for each site from cmdscale object and create a dataframe for plotting
pcoa.terr.bray.plotting <- as.data.frame(pcoa.terr.bray$points)
colnames(pcoa.terr.bray.plotting) <- c("axis_1", "axis_2")
pcoa.terr.bray.plotting$site <- rownames(pcoa.terr.bray.plotting)

# calculate the proportion of variance in the data which is explained by the first two PCoA axes
pcoa.terr.bray$eig[1]/(sum(pcoa.terr.bray$eig))
pcoa.terr.bray$eig[2]/(sum(pcoa.terr.bray$eig))

# create a PCoA plot
pcoa.terr.bray.plot <- ggplot(pcoa.terr.bray.plotting, aes(x = axis_1, y = axis_2, colour = site)) +
  geom_point(size = 3) +
  scale_colour_viridis_d(option = "magma", begin = 0.2, end = 0.8) +
  theme_bw() + 
  xlab("PCoA 1 (18.7%)") +
  ylab("PCoA 2 (16.7%)") +
  annotate(geom = 'text', label = 'Bray-Curtis', x = Inf, y = -Inf, hjust = 1.15, vjust = -1)

# repeat process with Jaccard dissimilarity matrix
pcoa.terr.jac <- cmdscale(lake.jac.terr, k = 2, eig = T)

pcoa.terr.jac.plotting <- as.data.frame(pcoa.terr.jac$points)
colnames(pcoa.terr.jac.plotting) <- c("axis_1", "axis_2")
pcoa.terr.jac.plotting$site <- rownames(pcoa.terr.jac.plotting)

pcoa.terr.jac$eig[1]/(sum(pcoa.terr.jac$eig))
pcoa.terr.jac$eig[2]/(sum(pcoa.terr.jac$eig))

pcoa.terr.jac.plot <- ggplot(pcoa.terr.jac.plotting, aes(x = axis_1, y = axis_2, colour = site)) +
  geom_point(size = 3) +
  scale_colour_viridis_d(option = "magma", begin = 0.2, end = 0.8) +
  theme_bw() + 
  xlab("PCoA 1 (16.2%)") +
  ylab("PCoA 2 (10.9%)") +
  annotate(geom = 'text', label = 'Jaccard', x = Inf, y = -Inf, hjust = 1.215, vjust = -1)

# extract plot legend
legend <- get_legend(pcoa.terr.jac.plot)

# plot Bray-Curtis PCoA and Jaccard PCoA side by side
pcoa.bray.jac.terr<-plot_grid(pcoa.terr.bray.plot + theme(legend.position = 'none'), pcoa.terr.jac.plot + theme(legend.position = 'none'), legend, ncol = 3, rel_widths = c(1,1,0.5))

##  PERMANOVA 

terr.permanova<-adonis2(lake.bray.terr ~ nhd_lat + julian.day + lake.area + IWS.stream.density + area.upstream.lakes + perc.agric.dev.land + max.depth,
       data = lake.env, permutations = 999)

##   ddRDA

dbRDA.terr.full <- capscale(lake.bray.terr ~ nhd_lat + julian.day+lake.area+IWS.stream.density+area.upstream.lakes+perc.agric.dev.land+max.depth,
                            lake.env.z)

vif.cca(dbRDA.terr.full)

anova(dbRDA.terr.full, by = "terms")
summary(dbRDA.terr.full)

smry.terr <- summary(dbRDA.terr.full)
scrs.terr <- scores(dbRDA.terr.full)
df.ddRAD.terr  <- data.frame(smry.terr$sites[,1:2]) # site scores for RDA1 and RDA2
df.ddRAD.terr$site <- rownames(df.ddRAD.terr)  #add site names
df.ddRAD.terr$site <-as.factor(df.ddRAD.terr$site)  #add site names

df.terr.merge<-cbind(df.ddRAD.terr,lake.env)
df.terr.merge<-df.terr.merge[c(1,2,3,6)]

df.ddRAD.terr2  <- data.frame(smry.terr$biplot[,1:2])  # mapping environmental variables

rda.terr.plot <- ggplot(df.terr.merge, aes(x=CAP1, y=CAP2, colour = nhd_lat)) + 
  geom_segment(data=df.ddRAD.terr2, aes(x=0, xend=CAP1, y=0, yend=CAP2), 
               color="grey50", arrow=arrow(length=unit(0.01,"npc"))) +
  geom_text(data=df.ddRAD.terr2, aes(x=CAP1,y=CAP2,label=rownames(df.ddRAD.terr2),
                                     hjust=0.5*(1-sign(CAP1)),vjust=0.5*(1-sign(CAP2))), 
            color="grey50", size=3) +
  geom_point(size = 3) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  xlim(-2, 2) +
  ylim(-2, 2) +
  xlab("RDA1 (36.2%)") + # this percentage comes from the CAP1 'importance of components:' proportion explained, which can be found in summary(dbRDA.mat) 
  ylab("RDA2 (20.1%)") + # this percentage comes from the CAP2 'importance of components:' proportion explained, which can be found in summary(dbRDA.mat) 
  coord_fixed() +
  theme_cowplot()+
  theme(legend.position = "none")
rda.terr.plot

#######################################
#         2.3 Wetland community       #  
#######################################

CM.filt.wet.lake<-read.csv("CM.filt.wet.lake.csv", header=TRUE)

df.wet<-CM.filt.wet.lake[,2:61]

lake.wet.mdf <- as.matrix.data.frame(df.wet)
rownames(lake.wet.mdf) <- lake.env$lake.name

lake.bray.wet <- vegdist(df.wet, method = "bray")
lake.bray.wet

lake.jac.wet <- vegdist(df.wet, method = "jaccard", binary = T)
lake.jac.wet

## Principal Coordinates Analysis (PCoA)

# calculate principal coordinates analysis (Bray-Curtis)
pcoa.wet.bray <- cmdscale(lake.bray.wet, k = 2, eig = T)

# extract axis positions for each site from cmdscale object and create a dataframe for plotting
pcoa.wet.bray.plotting <- as.data.frame(pcoa.wet.bray$points)
colnames(pcoa.wet.bray.plotting) <- c("axis_1", "axis_2")
pcoa.wet.bray.plotting$site <- rownames(pcoa.wet.bray.plotting)

# calculate the proportion of variance in the data which is explained by the first two PCoA axes
pcoa.wet.bray$eig[1]/(sum(pcoa.wet.bray$eig))
pcoa.wet.bray$eig[2]/(sum(pcoa.wet.bray$eig))

# create a PCoA plot
pcoa.wet.bray.plot <- ggplot(pcoa.wet.bray.plotting, aes(x = axis_1, y = axis_2, colour = site)) +
  geom_point(size = 3) +
  scale_colour_viridis_d(option = "magma", begin = 0.2, end = 0.8) +
  theme_bw() + 
  xlab("PCoA 1 (16.0%)") +
  ylab("PCoA 2 (11.8%)") +
  annotate(geom = 'text', label = 'Bray-Curtis', x = Inf, y = -Inf, hjust = 1.15, vjust = -1)

# repeat process with Jaccard dissimilarity matrix
pcoa.wet.jac <- cmdscale(lake.jac.wet, k = 2, eig = T)

pcoa.wet.jac.plotting <- as.data.frame(pcoa.wet.jac$points)
colnames(pcoa.wet.jac.plotting) <- c("axis_1", "axis_2")
pcoa.wet.jac.plotting$site <- rownames(pcoa.wet.jac.plotting)

pcoa.wet.jac$eig[1]/(sum(pcoa.wet.jac$eig))
pcoa.wet.jac$eig[2]/(sum(pcoa.wet.jac$eig))

pcoa.wet.jac.plot <- ggplot(pcoa.wet.jac.plotting, aes(x = axis_1, y = axis_2, colour = site)) +
  geom_point(size = 3) +
  scale_colour_viridis_d(option = "magma", begin = 0.2, end = 0.8) +
  theme_bw() + 
  xlab("PCoA 1 (21.4%)") +
  ylab("PCoA 2 (10.7%)") +
  annotate(geom = 'text', label = 'Jaccard', x = Inf, y = -Inf, hjust = 1.215, vjust = -1)

# extract plot legend
legend <- get_legend(pcoa.wet.jac.plot)

# plot Bray-Curtis PCoA and Jaccard PCoA side by side
pcoa.bray.jac.wet<-plot_grid(pcoa.wet.bray.plot + theme(legend.position = 'none'), pcoa.wet.jac.plot + theme(legend.position = 'none'), legend, ncol = 3, rel_widths = c(1,1,0.5))

## PERMANOVA

wet.permanova<-adonis2(lake.bray.wet ~ nhd_lat + julian.day + lake.area + IWS.stream.density + area.upstream.lakes + perc.agric.dev.land + max.depth,
       data = lake.env, permutations = 999)

## ddRDA

dbRDA.wet.full <- capscale(lake.bray.wet ~ nhd_lat + julian.day+lake.area+IWS.stream.density+area.upstream.lakes+perc.agric.dev.land+max.depth,
                           lake.env.z)

vif.cca(dbRDA.wet.full)

anova(dbRDA.wet.full, by = "terms")
summary(dbRDA.wet.full)

smry.wet <- summary(dbRDA.wet.full)
scrs.wet <- scores(dbRDA.wet.full)
df.ddRAD.wet  <- data.frame(smry.wet$sites[,1:2]) # site scores for RDA1 and RDA2
df.ddRAD.wet$site <- rownames(df.ddRAD.wet)  #add site names
df.ddRAD.wet$site <-as.factor(df.ddRAD.wet$site)  #add site names

df.wet.merge<-cbind(df.ddRAD.wet,lake.env)
df.wet.merge<-df.wet.merge[c(1,2,3,6)]

df.ddRAD.wet2  <- data.frame(smry.wet$biplot[,1:2])  # mapping environmental variables


rda.wet.plot<-ggplot(df.wet.merge, aes(x=CAP1, y=CAP2, colour = nhd_lat)) + 
  geom_segment(data=df.ddRAD.wet2, aes(x=0, xend=CAP1, y=0, yend=CAP2), 
               color="grey50", arrow=arrow(length=unit(0.01,"npc"))) +
  geom_text(data=df.ddRAD.wet2, aes(x=CAP1,y=CAP2,label=rownames(df.ddRAD.wet2),
                                    hjust=0.5*(1-sign(CAP1)),vjust=0.5*(1-sign(CAP2))), 
            color="grey50", size=3) +
  geom_point(size = 3) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  xlim(-2.5, 2) +
  ylim(-2, 1.5) +
  xlab("RDA1 (25.7%)") + # this percentage comes from the CAP1 'importance of components:' proportion explained, which can be found in summary(dbRDA.mat) 
  ylab("RDA2 (22.2%)") + # this percentage comes from the CAP2 'importance of components:' proportion explained, which can be found in summary(dbRDA.mat) 
  coord_fixed() +
  theme_cowplot()+
  theme(legend.position = "none")

rda.wet.plot

#######################################
#         2.4 Aquatic community       #  
#######################################

CM.filt.aqu.lake<-read.csv("CM.filt.aqu.lake.csv", header=TRUE)

df.aqu<-CM.filt.aqu.lake[,2:31]

lake.aqu.mdf <- as.matrix.data.frame(df.aqu)
rownames(lake.aqu.mdf) <- lake.env$lake.name

lake.bray.aqu <- vegdist(df.aqu, method = "bray")
lake.bray.aqu

lake.jac.aqu <- vegdist(df.aqu, method = "jaccard", binary = T)
lake.jac.aqu

## Principal Coordinates Analysis (PCoA)

# calculate principal coordinates analysis (Bray-Curtis)
pcoa.aqu.bray <- cmdscale(lake.bray.aqu, k = 2, eig = T)

# extract axis positions for each site from cmdscale object and create a dataframe for plotting
pcoa.aqu.bray.plotting <- as.data.frame(pcoa.aqu.bray$points)
colnames(pcoa.aqu.bray.plotting) <- c("axis_1", "axis_2")
pcoa.aqu.bray.plotting$site <- rownames(pcoa.aqu.bray.plotting)

# calculate the proportion of variance in the data which is explained by the first two PCoA axes
pcoa.aqu.bray$eig[1]/(sum(pcoa.aqu.bray$eig))
pcoa.aqu.bray$eig[2]/(sum(pcoa.aqu.bray$eig))

# create a PCoA plot
pcoa.aqu.bray.plot <- ggplot(pcoa.aqu.bray.plotting, aes(x = axis_1, y = axis_2, colour = site)) +
  geom_point(size = 3) +
  scale_colour_viridis_d(option = "magma", begin = 0.2, end = 0.8) +
  theme_bw() + 
  xlab("PCoA 1 (19.2%)") +
  ylab("PCoA 2 (18.1%)") +
  annotate(geom = 'text', label = 'Bray-Curtis', x = Inf, y = -Inf, hjust = 1.15, vjust = -1)

# repeat process with Jaccard dissimilarity matrix
pcoa.aqu.jac <- cmdscale(lake.jac.aqu, k = 2, eig = T)

pcoa.aqu.jac.plotting <- as.data.frame(pcoa.aqu.jac$points)
colnames(pcoa.aqu.jac.plotting) <- c("axis_1", "axis_2")
pcoa.aqu.jac.plotting$site <- rownames(pcoa.aqu.jac.plotting)

pcoa.aqu.jac$eig[1]/(sum(pcoa.aqu.jac$eig))
pcoa.aqu.jac$eig[2]/(sum(pcoa.aqu.jac$eig))

pcoa.aqu.jac.plot <- ggplot(pcoa.aqu.jac.plotting, aes(x = axis_1, y = axis_2, colour = site)) +
  geom_point(size = 3) +
  scale_colour_viridis_d(option = "magma", begin = 0.2, end = 0.8) +
  theme_bw() + 
  xlab("PCoA 1 (23.3%)") +
  ylab("PCoA 2 (17.8%)") +
  annotate(geom = 'text', label = 'Jaccard', x = Inf, y = -Inf, hjust = 1.215, vjust = -1)

# extract plot legend
legend <- get_legend(pcoa.aqu.jac.plot)

# plot Bray-Curtis PCoA and Jaccard PCoA side by side
pcoa.bray.jac.aqu<-plot_grid(pcoa.aqu.bray.plot + theme(legend.position = 'none'), pcoa.aqu.jac.plot + theme(legend.position = 'none'), legend, ncol = 3, rel_widths = c(1,1,0.5))

## PERMANOVA

aqu.permanova<-adonis2(lake.bray.aqu ~ nhd_lat + julian.day + lake.area + IWS.stream.density + area.upstream.lakes + perc.agric.dev.land + max.depth,
       data = lake.env, permutations = 999)

## ddRDA

dbRDA.aqu.full <- capscale(lake.bray.aqu ~ nhd_lat + julian.day+lake.area+IWS.stream.density+area.upstream.lakes+perc.agric.dev.land+max.depth,
                           lake.env.z)

vif.cca(dbRDA.aqu.full)

anova(dbRDA.aqu.full, by = "terms")
summary(dbRDA.aqu.full)

smry.aqu <- summary(dbRDA.aqu.full)
scrs.aqu <- scores(dbRDA.aqu.full)
df.ddRAD.aqu  <- data.frame(smry.aqu$sites[,1:2]) # site scores for RDA1 and RDA2
df.ddRAD.aqu$site <- rownames(df.ddRAD.aqu)  #add site names
df.ddRAD.aqu$site <-as.factor(df.ddRAD.aqu$site)  #add site names

df.aqu.merge<-cbind(df.ddRAD.aqu,lake.env)
df.aqu.merge<-df.aqu.merge[c(1,2,3,6)]

df.ddRAD.aqu2  <- data.frame(smry.aqu$biplot[,1:2])  # mapping environmental variables

rda.aqu.plot<-ggplot(df.aqu.merge, aes(x=CAP1, y=CAP2, colour = nhd_lat)) + 
  geom_segment(data=df.ddRAD.aqu2, aes(x=0, xend=CAP1, y=0, yend=CAP2), 
               color="grey50", arrow=arrow(length=unit(0.01,"npc"))) +
  geom_text(data=df.ddRAD.aqu2, aes(x=CAP1,y=CAP2,label=rownames(df.ddRAD.aqu2),
                                    hjust=0.5*(1-sign(CAP1)),vjust=0.5*(1-sign(CAP2))), 
            color="grey50", size=3) +
  geom_point(size = 3) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  xlim(-3, 3) +
  ylim(-2, 3) +
  xlab("RDA1 (34.0%)") + # this percentage comes from the CAP1 'importance of components:' proportion explained, which can be found in summary(dbRDA.mat) 
  ylab("RDA2 (21.1%)") + # this percentage comes from the CAP2 'importance of components:' proportion explained, which can be found in summary(dbRDA.mat) 
  coord_fixed() +
  theme_cowplot() 

rda.aqu.plot

combinedPlot <- plot_grid(rda.terr.plot,rda.wet.plot,rda.aqu.plot,
                           labels='auto', label_y= 1, label_x=0,align='hv',
                           label_size=20,ncol=3, rel_widths=1, rel_heights= 1)

png("dbRDA Figure.png",width = 4000, height = 1000, units = "px", res=300)
combinedPlot
dev.off()


###########################################
#             3. BETAPART                 #
###########################################

df.entire.betapart<-df.entire %>% mutate_if(is.numeric, ~1 * (. > 0))
entire.core<- betapart.core(df.entire.betapart)
entire.multi <- beta.multi(entire.core)

df.terr.betapart<-df.terr %>% mutate_if(is.numeric, ~1 * (. > 0))
terr.core<- betapart.core(df.terr.betapart)
terr.multi <- beta.multi(terr.core)

df.wet.betapart<-df.wet %>% mutate_if(is.numeric, ~1 * (. > 0))
wet.core<- betapart.core(df.wet.betapart)
wet.multi <- beta.multi(wet.core)

df.aqu.betapart<-df.aqu %>% mutate_if(is.numeric, ~1 * (. > 0))
aqu.core<- betapart.core(df.aqu.betapart)
aqu.multi <- beta.multi(aqu.core)



