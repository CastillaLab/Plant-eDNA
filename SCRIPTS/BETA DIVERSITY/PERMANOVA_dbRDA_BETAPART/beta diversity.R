############################################################
# Description:
#  - Filter samples based on sequencing depth
#  - Subset plant communities by habitat type
#  - Aggregate data at the lake level
#  - Perform PERMANOVA and dbRDA analyses
#  - Visualize ordination results
#  - Quantify beta diversity (turnover and nestedness)
############################################################

############################
#  LOAD REQUIRED PACKAGES  #
############################

library(ggplot2)  # Data visualization
library(cowplot)  # Plot formatting and combination
library(vegan)    # Community ecology analyses (PERMANOVA, dbRDA)
library(dplyr)    # Data manipulation
library(betapart) # Beta diversity partitioning


############################
#     DATA PREPARATION     #
############################

# Load plant community matrix
# Rows = samples; columns = species + metadata
CM<-read.csv("FINAL_plants.median.clean.CM_0diff_NO_gen.rev_blast_28Dec21.csv",header=TRUE)


######################################################
#      FILTER SAMPLES WITH LOW SEQUENCING DEPTH      #
######################################################

# Calculate total reads per sample (species columns assumed 7:396)
CM$total=rowSums(CM[,7:396])

# Inspect distribution of sequencing depth
hist(CM$total)
range(CM$total)

# Retain samples with >= 1000 reads
CM.filtered<-filter(CM,total >= 1000)

# Check filtered distribution
hist(CM.filtered$total)
range(CM.filtered$total)

####################################################
#          SUBSET SPECIES BY HABITAT TYPE          #
####################################################

# Load species metadata (species name + habitat classification)
spinfo<-read.csv("spinfo.csv", header=TRUE)

# Extract species lists by habitat
keep_spp_terr<-spinfo$species[spinfo$habitat == "terrestrial"]
keep_spp_wet<-spinfo$species[spinfo$habitat == "wetland"]
keep_spp_aqu<-spinfo$species[spinfo$habitat == "aquatic"]

# Subset community matrix for each habitat
CM.filt.terr<-CM.filtered[,c(1:3,which(colnames(CM.filtered) %in% keep_spp_terr))]
CM.filt.wet<-CM.filtered[,c(1:3,which(colnames(CM.filtered) %in% keep_spp_wet))]
CM.filt.aqu<-CM.filtered[,c(1:3,which(colnames(CM.filtered) %in% keep_spp_aqu))]



######################################
#       AGGREGATE DATA BY LAKE       #
######################################

# TERRESTRIAL
CM.filt.terr <- CM.filt.terr[,c(1,4:134)]
CM.filt.terr.lake <- CM.filt.terr %>%
  group_by(lake_name) %>%
  summarise(across(everything(), sum), .groups = "drop")

# WETLAND
CM.filt.wet <- CM.filt.wet[,c(1,4:65)]
CM.filt.wet.lake <- CM.filt.wet %>%
  group_by(lake_name) %>%
  summarise(across(everything(), sum), .groups = "drop")

# AQUATIC
CM.filt.aqu <- CM.filt.aqu[,c(1,4:31)]
CM.filt.aqu.lake <- CM.filt.aqu %>%
  group_by(lake_name) %>%
  summarise(across(everything(), sum), .groups = "drop")

###########################################
#    ENVIRONMENTAL DATA PREPARATION       #
###########################################

# Load environmental variables
lake.env<-read.csv("lakes traits.csv",header=TRUE)

# Copy for standardized version
lake.env.z <- lake.env

# Convert categorical variables to factors
lake.env.z$lake_name<-as.factor(lake.env.z$lake_name)
lake.env.z$year<-as.factor(lake.env.z$year)

# Standardize continuous variables (mean = 0, SD = 1)
lake.env.z$julian.day <- (lake.env.z$julian.day - mean(lake.env.z$julian.day))/sd(lake.env.z$julian.day)
lake.env.z$max.depth <- (lake.env.z$max.depth - mean(lake.env.z$max.depth))/sd(lake.env.z$max.depth)
lake.env.z$lake.area <- (lake.env.z$lake.area  - mean(lake.env.z$lake.area))/sd(lake.env.z$lake.area)
lake.env.z$IWS.stream.density <- (lake.env.z$IWS.stream.density - mean(lake.env.z$IWS.stream.density))/sd(lake.env.z$IWS.stream.density)
lake.env.z$area.upstream.lakes <- (lake.env.z$area.upstream.lakes - mean(lake.env.z$area.upstream.lakes))/sd(lake.env.z$area.upstream.lakes)
lake.env.z$perc.agric.dev.land <- (lake.env.z$perc.agric.dev.land - mean(lake.env.z$perc.agric.dev.land))/sd(lake.env.z$perc.agric.dev.land)
lake.env.z$nhd_lat <- (lake.env.z$nhd_lat - mean(lake.env.z$nhd_lat))/sd(lake.env.z$nhd_lat)


###################################################
#       TERRESTRIAL COMMUNITY ANALYSIS            #
###################################################

# Remove lake_name column for distance calculation
df.terr <- CM.filt.terr.lake[,2:ncol(CM.filt.terr.lake)]

# Compute Bray–Curtis dissimilarity
lake.bray.terr <- vegdist(df.terr, method = "bray")

##  PERMANOVA 
terr.permanova <- adonis2(
  lake.bray.terr ~ year + nhd_lat + julian.day + lake.area +
    IWS.stream.density + area.upstream.lakes +
    perc.agric.dev.land + max.depth,
  data = lake.env.z,
  permutations = 999,
  by = "terms")

terr.permanova


##  dbRDA
dbRDA.terr.full <- capscale(lake.bray.terr ~ year + nhd_lat + julian.day + 
                              lake.area + IWS.stream.density + area.upstream.lakes + 
                              perc.agric.dev.land + max.depth,
                            lake.env.z)

# Check multicollinearity
vif.cca(dbRDA.terr.full)

# Test significance of predictors
anova(dbRDA.terr.full, by = "terms")

# Model summary
summary(dbRDA.terr.full)

##############################
# PREPARE DATA FOR PLOTTING  #
##############################

smry.terr <- summary(dbRDA.terr.full)


# Site scores (ordination coordinates)
df.ddRAD.terr  <- data.frame(smry.terr$sites[,1:2])
df.ddRAD.terr$site <- rownames(df.ddRAD.terr)  #add site names

# Merge with environmental data
df.terr.merge<-cbind(df.ddRAD.terr,lake.env)
df.terr.merge<-df.terr.merge[c(1,2,3,15)]

# Environmental vectors
df.ddRAD.terr2  <- data.frame(smry.terr$biplot[,1:2])  # mapping environmental variables


#####################
#     PLOT dbRDA    #
#####################
rda.terr.plot <- ggplot(df.terr.merge, aes(x = CAP1, y = CAP2, colour = nhd_lat)) +
  geom_segment(data = df.ddRAD.terr2,
               aes(x = 0, xend = CAP1, y = 0, yend = CAP2),
               color = "grey50",
               arrow = arrow(length = unit(0.01, "npc"))) +
  geom_text(data = df.ddRAD.terr2,
            aes(label = rownames(df.ddRAD.terr2)),
            color = "grey50", size = 3) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  theme_cowplot() +
  theme(legend.position = "none")

rda.terr.plot


#######################################
#     WETLAND COMMUNITY ANALYSIS     #
#######################################

df.wet <- CM.filt.wet.lake[,2:ncol(CM.filt.wet.lake)]
lake.bray.wet <- vegdist(df.wet, method = "bray")



# PERMANOVA
wet.permanova <- adonis2(lake.bray.wet ~ year + nhd_lat + julian.day + lake.area +
                           IWS.stream.density + area.upstream.lakes +
                           perc.agric.dev.land + max.depth,
                         data = lake.env.z, permutations = 999, by = "terms")

wet.permanova

##  dbRDA
dbRDA.wet.full <- capscale(lake.bray.wet ~ year + nhd_lat + julian.day +
                             lake.area + IWS.stream.density + area.upstream.lakes +
                             perc.agric.dev.land + max.depth,
                           data = lake.env.z)

# Check multicollinearity
vif.cca(dbRDA.wet.full)

# Test significance of predictors
anova(dbRDA.wet.full, by = "terms")

# Model summary
summary(dbRDA.wet.full)


##############################
# PREPARE DATA FOR PLOTTING  #
##############################
smry.wet <- summary(dbRDA.wet.full)

# Site scores (ordination coordinates)
df.ddRAD.wet  <- data.frame(smry.wet$sites[,1:2]) 
df.ddRAD.wet$site <- rownames(df.ddRAD.wet) 

# Merge with environmental data
df.wet.merge<-cbind(df.ddRAD.wet,lake.env)
df.wet.merge<-df.wet.merge[c(1,2,3,15)]

# Environmental vectors
df.ddRAD.wet2  <- data.frame(smry.wet$biplot[,1:2])  # mapping environmental variables


#####################
#     PLOT dbRDA    #
#####################

rda.wet.plot<-ggplot(df.wet.merge, aes(x=CAP1, y=CAP2, colour = nhd_lat)) + 
  geom_segment(data=df.ddRAD.wet2, 
               aes(x=0, xend=CAP1, y=0, yend=CAP2), 
               color="grey50", 
               arrow=arrow(length=unit(0.01,"npc"))) +
  geom_text(data=df.ddRAD.wet2, 
            aes(label=rownames(df.ddRAD.wet2)),
                color="grey50", size=3) +
  geom_point(size = 3) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  theme_cowplot()+
  theme(legend.position = "none")

rda.wet.plot


##########################################
#       AQUATIC COMMUNITY ANALYSIS       #
##########################################

df.aqu <- CM.filt.aqu.lake[,2:ncol(CM.filt.aqu.lake)]
lake.bray.aqu <- vegdist(df.aqu, method = "bray")



## PERMANOVA
aqu.permanova<-adonis2(lake.bray.aqu ~ year + nhd_lat + julian.day + 
                         lake.area + IWS.stream.density + 
                         area.upstream.lakes + perc.agric.dev.land + 
                         max.depth, data = lake.env.z, permutations = 999,
                       by="terms")
aqu.permanova


## ddRDA

dbRDA.aqu.full <- capscale(lake.bray.aqu ~ year + nhd_lat + julian.day + 
                             lake.area + IWS.stream.density + area.upstream.lakes + 
                             perc.agric.dev.land + max.depth,lake.env.z)

# Check multicollinearity
vif.cca(dbRDA.aqu.full)

# Test significance of predictors
anova(dbRDA.aqu.full, by = "terms")

# Model summary
summary(dbRDA.aqu.full)


##############################
# PREPARE DATA FOR PLOTTING  #
##############################
smry.aqu <- summary(dbRDA.aqu.full)

# Site scores (ordination coordinates)
df.ddRAD.aqu  <- data.frame(smry.aqu$sites[,1:2]) 
df.ddRAD.aqu$site <- rownames(df.ddRAD.aqu)  

# Merge with environmental data
df.aqu.merge<-cbind(df.ddRAD.aqu,lake.env)
df.aqu.merge<-df.aqu.merge[c(1,2,3,15)]

# Environmental vectors
df.ddRAD.aqu2  <- data.frame(smry.aqu$biplot[,1:2])  # mapping environmental variables


#####################
#     PLOT dbRDA    #
#####################
rda.aqu.plot<-ggplot(df.aqu.merge, aes(x=CAP1, y=CAP2, colour = nhd_lat)) + 
  geom_segment(data=df.ddRAD.aqu2, 
               aes(x=0, xend=CAP1, y=0, yend=CAP2), 
               color="grey50", 
               arrow=arrow(length=unit(0.01,"npc"))) +
  geom_text(data=df.ddRAD.aqu2, 
            aes(label=rownames(df.ddRAD.aqu2)),
            color="grey50", size=3) +
  geom_point(size = 3) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  theme_cowplot() 

rda.aqu.plot


combinedPlot <- plot_grid(rda.terr.plot,rda.wet.plot,rda.aqu.plot,
                           labels='auto', label_y= 1, label_x=0,align='hv',
                           label_size=20,ncol=3, rel_widths=1, rel_heights= 1)

png("dbRDA Figure.png",width = 4000, height = 1000, units = "px", res=300)
combinedPlot
dev.off()


###########################################
#                BETAPART                 #
###########################################

# Purpose:
# This section computes beta diversity partitioning
# (turnover and nestedness components) using the
# betapart package for three habitat types:
# terrestrial, wetland, and aquatic communities.
#
# Steps:
# 1. Convert abundance/continuous data to presence–absence (0/1)
# 2. Build betapart core objects
# 3. Compute multiple-site beta diversity indices
###########################################

# ---------------------------
# TERRESTRIAL COMMUNITIES
# ---------------------------

# Convert all numeric columns to presence–absence (1 = present, 0 = absent)
df.terr.betapart<-df.terr %>% mutate_if(is.numeric, ~1 * (. > 0))

# Create betapart core object
terr.core<- betapart.core(df.terr.betapart)

# Calculate multiple-site beta diversity
terr.multi <- beta.multi(terr.core)

# ---------------------------
# WETLAND COMMUNITIES
# ---------------------------

# Convert to presence–absence matrix
df.wet.betapart<-df.wet %>% mutate_if(is.numeric, ~1 * (. > 0))

# Build betapart core object
wet.core<- betapart.core(df.wet.betapart)

# Compute beta diversity components
wet.multi <- beta.multi(wet.core)

# ---------------------------
# AQUATIC COMMUNITIES
# ---------------------------

# Convert to presence–absence matrix
df.aqu.betapart<-df.aqu %>% mutate_if(is.numeric, ~1 * (. > 0))

# Build betapart core object
aqu.core<- betapart.core(df.aqu.betapart)

# Compute beta diversity components
aqu.multi <- beta.multi(aqu.core)



