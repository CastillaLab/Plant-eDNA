############################################################
# ANALYSIS OF PLANT COMMUNITY COMPOSITION ACROSS LAKES
# - Filters sequencing data
# - Subsets species by habitat type
# - Aggregates reads per lake
# - Computes proportional composition
# - Fits beta regression models (glmmTMB)
# - Generates prediction plots
############################################################


########################
# Load packages
########################
# Data wrangling and plotting
library(dplyr)
library(tidyverse)
library(ggplot2)
library(cowplot)

# Modeling framework (beta regression)
library(glmmTMB)

# Modeling framework (beta regression)
library(lmerTest) 
library(lmtest)
library(MuMIn)

# Collinearity diagnostics
library(car)

# Predicted values for plotting
library(ggeffects)

########################
# Load community matrix
########################
# Rows = samples (lakes), columns = taxa + metadata
CM<-read.csv("FINAL_plants.median.clean.CM_0diff_NO_gen.rev_blast_28Dec21.csv",header=TRUE)


######################################################
# FILTER SAMPLES WITH LOW SEQUENCE DEPTH
######################################################

# Compute total reads per sample (sum across species columns)
CM$total=rowSums(CM[,7:396])

# Inspect distribution of sequencing depth
hist(CM$total)
range(CM$total)

# Keep only samples with sufficient sequencing depth (>= 1000 reads)
CM.filtered<-filter(CM,total >= 1000)

# Verify filtering outcome
hist(CM.filtered$total)
range(CM.filtered$total)

####################################################
# SUBSET SPECIES BY HABITAT TYPE
####################################################

# Load species metadata (habitat classification)
spinfo<-read.csv("spinfo.csv", header=TRUE)

# -------------------------
# Terrestrial species
# -------------------------
keep_spp_terr<-spinfo$species[spinfo$habitat == "terrestrial"]
CM.filt.terr<-CM.filtered[,c(1:3,which(colnames(CM.filtered) %in% keep_spp_terr))]

# -------------------------
# Wetland species
# -------------------------
keep_spp_wet<-spinfo$species[spinfo$habitat == "wetland"]
CM.filt.wet<-CM.filtered[,c(1:3,which(colnames(CM.filtered) %in% keep_spp_wet))]

# -------------------------
# Aquatic species
# -------------------------
keep_spp_aqu<-spinfo$species[spinfo$habitat == "aquatic"]
CM.filt.aqu<-CM.filtered[,c(1:3,which(colnames(CM.filtered) %in% keep_spp_aqu))]


####################################################
# AGGREGATE READS PER LAKE
####################################################

# Reduce each habitat matrix to species columns only
CM.filt.terr<-CM.filt.terr[,c(1,4:134)]
CM.filt.wet<-CM.filt.wet[,c(1,4:65)]
CM.filt.aqu<-CM.filt.aqu[,c(1,4:31)]

# Sum reads per lake for each species (aggregation step)
CM.filt.terr.lake<-CM.filt.terr %>% 
  group_by(lake_name) %>% summarise(across(everything(),list(sum)))

CM.filt.wet.lake<-CM.filt.wet %>% 
  group_by(lake_name) %>% summarise(across(everything(),list(sum)))

CM.filt.aqu.lake<-CM.filt.aqu %>% 
  group_by(lake_name) %>% summarise(across(everything(),list(sum)))

####################################################
# TOTAL READS PER HABITAT AND LAKE
####################################################

# Collapse species columns into total reads per habitat

CM.filt.terr.lake$terrestrial.reads=rowSums(CM.filt.terr.lake[,2:132])
CM.filt.wet.lake$wetland.reads=rowSums(CM.filt.wet.lake[,2:63])
CM.filt.aqu.lake$aquatic.reads=rowSums(CM.filt.aqu.lake[,2:29])

####################################################
# COMBINE HABITATS INTO SINGLE DATA FRAME
####################################################

df_list<-list(CM.filt.terr.lake,CM.filt.wet.lake,CM.filt.aqu.lake)

# Merge all habitat-specific summaries by lake
df.prop<-df_list %>% 
  reduce(full_join, by='lake_name')

# Keep only total habitat reads
df.prop<-df.prop %>% 
  dplyr::select(lake_name,terrestrial.reads,wetland.reads,aquatic.reads)

# Compute total reads per lake
df.prop$total.reads<-df.prop$terrestrial.reads +
  df.prop$wetland.reads + 
  df.prop$aquatic.reads

# Convert to proportional composition
df.prop$prop.terrestrial<-df.prop$terrestrial.reads/df.prop$total.reads
df.prop$prop.wetland<-df.prop$wetland.reads/df.prop$total.reads
df.prop$prop.aquatic<-df.prop$aquatic.reads/df.prop$total.reads

####################################################
# MERGE WITH LAKE AND LANDSCAPE TRAITS
####################################################

lake.traits<-read.csv("lakes traits.csv", header=TRUE)
merged.df<-merge(df.prop,lake.traits,by="lake_name")


####################################################
# CREATE DERIVED PREDICTOR VARIABLES
####################################################

merged.df$forest.500 <- merged.df$perc_forest_2011_buffer500
merged.df$developed.500 <- merged.df$perc_dev_land_2011_buffer500
merged.df$agric.500 <- merged.df$perc_agric_land_2011_buffer500
merged.df$log.lake.area <- log(merged.df$lake.area)

####################################################
# MULTICOLLINEARITY DIAGNOSTIC (VIF)
####################################################

vif_model_terr <- lm(prop.terrestrial ~ log.lake.area + max.depth + forest.500 + agric.500 + developed.500,
                     data = merged.df)
vif(vif_model_terr)


####################################################
# MODEL SELECTION (TERRESTRIAL)
####################################################

# Candidate beta regression models for terrestrial proportion
m1a <- glmmTMB(prop.terrestrial ~ log.lake.area, family=beta_family(), data=merged.df)
m2a <- glmmTMB(prop.terrestrial ~ max.depth, family=beta_family(), data=merged.df)
m3a <- glmmTMB(prop.terrestrial ~ forest.500, family=beta_family(), data=merged.df)
m4a <- glmmTMB(prop.terrestrial ~ developed.500, family=beta_family(), data=merged.df)
m5a <- glmmTMB(prop.terrestrial ~ agric.500, family=beta_family(), data=merged.df)
m6a <- glmmTMB(prop.terrestrial ~ log.lake.area + max.depth + forest.500, family=beta_family(), data=merged.df)
m7a <- glmmTMB(prop.terrestrial ~ agric.500 + developed.500, family=beta_family(), data=merged.df)
m8a <- glmmTMB(prop.terrestrial ~ log.lake.area + max.depth + forest.500 + agric.500 + developed.500, 
           family=beta_family(), data=merged.df)
m9a <- glmmTMB(prop.terrestrial ~ 1, family=beta_family(), data=merged.df)

# Model comparison
output.terr <- model.sel(m1a,m2a,m3a,m4a,m5a,m6a,m7a,m8a,m9a)
output.terr

# Likelihood ratio tests for selected comparisons
lrtest(m3a,m9a)
lrtest(m4a,m9a)


####################################################
# MODEL SELECTION (WETLAND)
####################################################

# Full model used for diagnostics

full.wet <- glmmTMB(prop.wetland ~ log.lake.area + max.depth + 
                  forest.500 + agric.500 + developed.500,
                family = beta_family(), data = merged.df)

# VIF diagnostics
vif_model_wet <- lm(prop.wetland ~ log.lake.area + max.depth + forest.500 + agric.500 + developed.500,
                     data = merged.df)
vif(vif_model_wet)

# Candidate models
m1b <- glmmTMB(prop.wetland ~ log.lake.area, family=beta_family(), data=merged.df)
m2b <- glmmTMB(prop.wetland ~ max.depth, family=beta_family(), data=merged.df)
m3b <- glmmTMB(prop.wetland ~ forest.500, family=beta_family(), data=merged.df)
m4b <- glmmTMB(prop.wetland ~ developed.500, family=beta_family(), data=merged.df)
m5b <- glmmTMB(prop.wetland ~ agric.500, family=beta_family(), data=merged.df)
m6b <- glmmTMB(prop.wetland ~ log.lake.area + max.depth + forest.500, family=beta_family(), data=merged.df)
m7b <- glmmTMB(prop.wetland ~ agric.500 + developed.500, family=beta_family(), data=merged.df)
m8b <- glmmTMB(prop.wetland ~ log.lake.area + max.depth + forest.500 + agric.500 + developed.500, 
           family=beta_family(), data=merged.df)
m9b <- glmmTMB(prop.wetland ~ 1, family=beta_family(), data=merged.df)

output.wet <- model.sel(m1b,m2b,m3b,m4b,m5b,m6b,m7b,m8b,m9b)
output.wet

lrtest(m3b,m9b)


####################################################
# MODEL SELECTION (AQUATIC)
####################################################

# Full model + VIF diagnostics

full.aqu <- glmmTMB(prop.aquatic ~ log.lake.area + max.depth + 
                  forest.500 + agric.500 + developed.500,
                family = beta_family(), data = merged.df)

vif_model_aqu <- lm(prop.aquatic ~ log.lake.area + max.depth + forest.500 + agric.500 + developed.500,
                    data = merged.df)
vif(vif_model_aqu)

# Candidate models
m1c <- glmmTMB(prop.aquatic ~ log.lake.area, family=beta_family(), data=merged.df)
m2c <- glmmTMB(prop.aquatic ~ max.depth, family=beta_family(), data=merged.df)
m3c <- glmmTMB(prop.aquatic ~ forest.500, family=beta_family(), data=merged.df)
m4c <- glmmTMB(prop.aquatic ~ developed.500, family=beta_family(), data=merged.df)
m5c <- glmmTMB(prop.aquatic ~ agric.500, family=beta_family(), data=merged.df)
m6c <- glmmTMB(prop.aquatic ~ log.lake.area + max.depth + forest.500, family=beta_family(), data=merged.df)
m7c <- glmmTMB(prop.aquatic ~ agric.500 + developed.500, family=beta_family(), data=merged.df)
m8c <- glmmTMB(prop.aquatic ~ log.lake.area + max.depth + forest.500 + agric.500 + developed.500, 
           family=beta_family(), data=merged.df)
m9c <- glmmTMB(prop.aquatic ~ 1, family=beta_family(), data=merged.df)

output.aqu <- model.sel(m1c,m2c,m3c,m4c,m5c,m6c,m7c,m8c,m9c)
output.aqu

lrtest(m4c,m9c)
lrtest(m3c,m9c)

####################################################
# PREDICTIONS FOR PLOTTING
####################################################

# Generate fitted values for observed data
merged.df$pred_m4a <- predict(m4a, type = "response")
merged.df$pred_m3a <- predict(m3a, type = "response")
merged.df$pred_m3b <- predict(m3b, type = "response")
merged.df$pred_m4c <- predict(m4c, type = "response")
merged.df$pred_m3c <- predict(m3c, type = "response")

####################################################
# COMMON PLOT THEME
####################################################

common_theme <- theme(
  legend.title = element_blank(),
  panel.background = element_rect(fill = "white"),
  panel.border = element_rect(linetype = "solid", fill = NA),
  axis.text = element_text(size = 14),
  axis.title.y = element_text(margin = margin(r = 20), size = 16),
  axis.title.x = element_text(margin = margin(t = 20), size = 16),
  plot.margin = margin(10, 10, 10, 30)
)

####################################################
# PLOTTING SECTION (MULTIPLE LAND USE EFFECTS)
####################################################

# Each plot shows observed values, predicted values, and model fit from ggpredict

#############################
#    Terrestrial species    #
#############################

# Developed land
m4a.pred <- ggpredict(m4a, terms = "developed.500")

p1<-ggplot() +
  geom_line(data = m4a.pred, aes(x = x, y = predicted),
            linewidth = 1.5, color = "darkgray") +
  geom_ribbon(data = m4a.pred,
              aes(x = x, ymin = conf.low, ymax = conf.high),
              alpha = 0.2) +
  geom_point(data = merged.df,
             aes(x = developed.500, y = prop.terrestrial),
             size = 5, shape = 21, stroke = 0.8) +
  geom_point(data = merged.df,
             aes(x = developed.500, y = pred_m4a),
             fill = "coral", color = "black",
             shape = 21, size = 5) +
  ylim(0,1) +
  common_theme +
  ylab("Proportion of\n terrestrial species") +
  xlab("% Developed land (500 m)")

# Forest cover
m3a.pred <- ggpredict(m3a, terms = "forest.500")

p2 <- ggplot() +
  geom_line(data = m3a.pred, aes(x = x, y = predicted),
            linewidth = 1.5, color = "darkgray") +
  geom_ribbon(data = m3a.pred,
              aes(x = x, ymin = conf.low, ymax = conf.high),
              alpha = 0.2) +
  geom_point(data = merged.df,
             aes(x = forest.500, y = prop.terrestrial),
             size = 5, shape = 21, stroke = 0.8) +
  geom_point(data = merged.df,
             aes(x = forest.500, y = pred_m3a),
             fill = "coral", color = "black",
             shape = 21, size = 5) +
  ylim(0,1) +
  common_theme +
  ylab("Proportion of\n terrestrial species") +
  xlab("% Forest cover (500 m)")

#########################
#    Wetland species    #
#########################

# Forest cover
m3b.pred <- ggpredict(m3b, terms = "forest.500")

p3 <- ggplot() +
  geom_line(data = m3b.pred, aes(x = x, y = predicted),
            linewidth = 1.5, color = "darkgray") +
  geom_ribbon(data = m3b.pred,
              aes(x = x, ymin = conf.low, ymax = conf.high),
              alpha = 0.2) +
  geom_point(data = merged.df,
             aes(x = forest.500, y = prop.wetland),
             size = 5, shape = 21, stroke = 0.8) +
  geom_point(data = merged.df,
             aes(x = forest.500, y = pred_m3b),
             fill = "coral", color = "black",
             shape = 21, size = 5) +
  ylim(0,1) +
  common_theme +
  ylab("Proportion of\n wetland species") +
  xlab("% Forest cover (500 m)")

#########################
#    Aquatic species    #
#########################

# Developed land
m4c.pred <- ggpredict(m4c, terms = "developed.500")

p4 <- ggplot() +
  geom_line(data = m4c.pred, aes(x = x, y = predicted),
            linewidth = 1.5, color = "darkgray") +
  geom_ribbon(data = m4c.pred,
              aes(x = x, ymin = conf.low, ymax = conf.high),
              alpha = 0.2) +
  geom_point(data = merged.df,
             aes(x = developed.500, y = prop.aquatic),
             size = 5, shape = 21, stroke = 0.8) +
  geom_point(data = merged.df,
             aes(x = developed.500, y = pred_m4c),
             fill = "coral", color = "black",
             shape = 21, size = 5) +
  ylim(0,1) +
  common_theme +
  ylab("Proportion of\n aquatic species") +
  xlab("% Developed land (500 m)")

# Forest cover
m3c.pred <- ggpredict(m3c, terms = "forest.500")

p5<-ggplot() +
  geom_line(data = m3c.pred, aes(x = x, y = predicted),
            linewidth = 1.5, color = "darkgray") +
  geom_ribbon(data = m3c.pred,
              aes(x = x, ymin = conf.low, ymax = conf.high),
              alpha = 0.2) +
  geom_point(data = merged.df,
             aes(x = forest.500, y = prop.aquatic),
             size = 5, shape = 21, stroke = 0.8) +
  geom_point(data = merged.df,
             aes(x = forest.500, y = pred_m3c),
             fill = "coral", color = "black",
             shape = 21, size = 5) +
  ylim(0,1) +
  common_theme +
  ylab("Proportion of\n aquatic species") +
  xlab("% Forest (500 m)")


####################################################
# COMBINE PLOTS
####################################################

combinedPlot1 <- plot_grid(
  p1, p2, p3, p4, p5,
  labels = "auto",
  label_y = 1,
  label_x = 0.05,
  align = "hv",
  label_size = 20,
  ncol = 3,
  nrow = 2
)

####################################################
# EXPORT FIGURE
####################################################

png("Proportions regression models.png",
    width = 4500, height = 2200,
    units = "px", res = 300)

combinedPlot1
dev.off()

