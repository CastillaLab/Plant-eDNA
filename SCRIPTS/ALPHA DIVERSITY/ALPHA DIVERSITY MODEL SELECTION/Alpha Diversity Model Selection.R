############################################################
# Description:
#  - Load plant community matrix and species metadata
#  - Filter low-read samples (<1000 reads)
#  - Subset species by habitat (terrestrial, wetland, aquatic)
#  - Aggregate data at the lake level
#  - Compute alpha diversity metrics (Shannon and richness)
#  - Fit GLM / negative binomial models
#  - Perform model selection and likelihood ratio tests
#  - Generate model predictions and visualization
############################################################


############################
## LOAD REQUIRED PACKAGES ##
############################

library(tidyverse)    # data manipulation and plotting
library(vegan)        # diversity metrics (Shannon, richness)
library(MuMIn)        # model selection (AIC-based)
library(MASS)         # negative binomial models (glm.nb)
library(lmtest)       # likelihood ratio tests
library(ggeffects)    # model predictions
library(car)          # multicollinearity (VIF)
library(cowplot)      # multi-panel plots



############################
##   1. DATA PREPARATION  ##  
############################

# Load the plant community matrix
CM<-read.csv("FINAL_plants.median.clean.CM_0diff_NO_gen.rev_blast_28Dec21.csv",header=TRUE)


######################################################
# FILTER SAMPLES BASED ON SEQUENCING DEPTH
######################################################

# Compute total reads per sample (species columns assumed 7:396)
CM$total=rowSums(CM[,7:396])

# Inspect distribution of sequencing depth
hist(CM$total)
range(CM$total)

# Retain samples with sufficient reads (>= 1000)
CM.filtered<-filter(CM,total >= 1000)

# Re-check after filtering
hist(CM.filtered$total)
range(CM.filtered$total)

#####################################################
# SUBSET SPECIES BY HABITAT TYPE
#####################################################

# Load species metadata (includes habitat classification)
spinfo<-read.csv("spinfo.csv", header=TRUE)

#######################
# TERRESTRIAL SPECIES #
#######################

# --- Terrestrial species ---
keep_spp_terr<-spinfo$species[spinfo$habitat == "terrestrial"]
CM.filt.terr<-CM.filtered[,c(1:3,which(colnames(CM.filtered) %in% keep_spp_terr))]

# --- Wetland species ---
keep_spp_wet<-spinfo$species[spinfo$habitat == "wetland"]
CM.filt.wet<-CM.filtered[,c(1:3,which(colnames(CM.filtered) %in% keep_spp_wet))]

# --- Aquatic species ---
keep_spp_aqu<-spinfo$species[spinfo$habitat == "aquatic"]
CM.filt.aqu<-CM.filtered[,c(1:3,which(colnames(CM.filtered) %in% keep_spp_aqu))]

######################################
## AGGREGATE DATA AT LAKE LEVEL     ##
######################################

# Sum species abundances within each lake

# Terrestrial
CM.filt.terr<-CM.filt.terr[,c(1,4:134)]
CM.filt.terr.lake<-CM.filt.terr %>% 
  group_by(lake_name) %>% 
  summarise(across(everything(),list(sum)))


# Wetland
CM.filt.wet<-CM.filt.wet[,c(1,4:65)]
CM.filt.wet.lake<-CM.filt.wet %>% 
  group_by(lake_name) %>% 
  summarise(across(everything(),list(sum)))


# Aquatic
CM.filt.aqu<-CM.filt.aqu[,c(1,4:31)]
CM.filt.aqu.lake<-CM.filt.aqu %>% 
  group_by(lake_name) %>% 
  summarise(across(everything(),list(sum)))


#####################################
## 2. ALPHA DIVERSITY CALCULATION  ##
#####################################

# Load lake-level environmental predictors
lake_traits<-read.csv("lakes traits.csv", header = T)

################################
## 2.1 TERRESTRIAL COMMUNITY  ##
################################

# Extract species matrix
df.terr<-CM.filt.terr.lake[,2:132]

# Compute diversity metrics

H.terr <- diversity(df.terr)            # Shannon diversity
richness.terr <- specnumber(df.terr)    # Species richness


# Combine with lake metadata
terr.mod.df <- merge(lake_traits,
                     cbind(lake_name = CM.filt.terr.lake$lake_name,
                           H.terr, richness.terr),
                     by = "lake_name")


################################
## 2.2 WETLAND COMMUNITY      ##
################################

df.wet<-CM.filt.wet.lake[,2:63]


H.wet <- diversity(df.wet)
richness.wet <- specnumber(df.wet)  

wet.mod.df <- merge(lake_traits,
                    cbind(lake_name = CM.filt.wet.lake$lake_name,
                          H.wet, richness.wet),
                    by = "lake_name")

################################
## 2.3 AQUATIC COMMUNITY      ##
################################

df.aqu<-CM.filt.aqu.lake[,2:29]


H.aqu <- diversity(df.aqu)
richness.aqu <- specnumber(df.aqu)  

aqu.mod.df <- merge(lake_traits,
                    cbind(lake_name = CM.filt.aqu.lake$lake_name,
                          H.aqu, richness.aqu),
                    by = "lake_name")


#####################################
#      3. MODELING APPROACH         #
#####################################

# General strategy:
#  - Transform predictors (log)
#  - Fit Gaussian models for Shannon diversity
#  - Fit negative binomial models for richness
#  - Compare models using AIC (MuMIn::model.sel)
#  - Evaluate significance with LR tests


##################################
# 3.1 TERRESTRIAL MODELS         #
##################################

# Transform predictors
terr.mod.df$log.lake.area<-log(terr.mod.df$lake.area)
terr.mod.df$log.area.upstream.lakes<-log(terr.mod.df$area.upstream.lakes+1)

# Ensure correct data types
terr.mod.df$H.terr<-as.numeric(terr.mod.df$H.terr)
terr.mod.df$richness.terr<-as.integer(terr.mod.df$richness.terr)
terr.mod.df$year<-as.factor(terr.mod.df$year)

# Full models (used to assess multicollinearity)
full.terr.H<-glm(H.terr ~ log.lake.area + max.depth + IWS.stream.density + log.area.upstream.lakes + 
              perc.agric.dev.land + julian.day, family = gaussian, data = terr.mod.df)

vif(full.terr.H)

full.terr.rich<-glm.nb(richness.terr~ log.lake.area+max.depth+IWS.stream.density+log.area.upstream.lakes+
                 perc.agric.dev.land+julian.day, data = terr.mod.df)
vif(full.terr.rich)


# Candidate models (single predictors + full + null)
# Used for model selection via AIC
# (same structure repeated for wetland and aquatic)

# Candidate models: Shannon Diversity - Terrestrial species
m1a<-glm(H.terr~year, family=gaussian, data=terr.mod.df)
m2a<-glm(H.terr~log.lake.area, family=gaussian, data=terr.mod.df)
m3a<-glm(H.terr~max.depth, family=gaussian, data=terr.mod.df)
m4a<-glm(H.terr~IWS.stream.density, family=gaussian, data=terr.mod.df)
m5a<-glm(H.terr~log.area.upstream.lakes, family=gaussian, data=terr.mod.df)
m6a<-glm(H.terr~perc.agric.dev.land, family=gaussian, data=terr.mod.df)
m7a<-glm(H.terr~julian.day, family=gaussian, data=terr.mod.df)
m8a<-glm(H.terr~ log.lake.area+max.depth+IWS.stream.density+log.area.upstream.lakes+perc.agric.dev.land+julian.day+year, family = gaussian, data = terr.mod.df)
m9a<-glm(H.terr~1,family=gaussian, data=terr.mod.df)

# Model selection
output.terr.H<-model.sel(m1a,m2a,m3a,m4a,m5a,m6a,m7a,m8a,m9a)

# Likelihood ratio test
lrtest(m7a,m9a)


# Candidate models: Species richness - Terrestrial species
m1b<-glm.nb(richness.terr~year, data=terr.mod.df)
m2b<-glm.nb(richness.terr~log.lake.area, data=terr.mod.df)
m3b<-glm.nb(richness.terr~max.depth, data=terr.mod.df)
m4b<-glm.nb(richness.terr~IWS.stream.density, data=terr.mod.df)
m5b<-glm.nb(richness.terr~log.area.upstream.lakes, data=terr.mod.df)
m6b<-glm.nb(richness.terr~perc.agric.dev.land, data=terr.mod.df)
m7b<-glm.nb(richness.terr~julian.day, data=terr.mod.df)
m8b<-glm.nb(richness.terr~ log.lake.area+max.depth+IWS.stream.density+log.area.upstream.lakes+perc.agric.dev.land+julian.day+year, data = terr.mod.df)
m9b<-glm.nb(richness.terr~1, data=terr.mod.df)

# Model selection
output.terr.rich<-model.sel(m1b,m2b,m3b,m4b,m5b,m6b,m7b,m8b,m9b)

# Likelihood ratio test
lrtest(m1b,m9b)
lrtest(m7b,m9b)

##################################
#     3.2 WETLAND MODELS         #
##################################

# Transform predictors
wet.mod.df$log.lake.area<-log(wet.mod.df$lake.area)
wet.mod.df$log.area.upstream.lakes<-log(wet.mod.df$area.upstream.lakes+1)

# Ensure correct data types
wet.mod.df$H.wet<-as.numeric(wet.mod.df$H.wet)
wet.mod.df$richness.wet<-as.integer(wet.mod.df$richness.wet)
wet.mod.df$year<-as.factor(wet.mod.df$year)

# Full models (used to assess multicollinearity)
full.wet.H<-glm(H.wet ~ log.lake.area + max.depth + IWS.stream.density + log.area.upstream.lakes + 
              perc.agric.dev.land + julian.day, family = gaussian, data = wet.mod.df)

vif(full.wet.H)

full.wet.rich<-glm.nb(richness.wet~ log.lake.area+max.depth+IWS.stream.density+log.area.upstream.lakes+
                 perc.agric.dev.land+julian.day, data = wet.mod.df)
vif(full.wet.rich)


# Candidate models: Shannon Diversity - Wetland species

m1c<-glm(H.wet~year, family=gaussian, data=wet.mod.df)
m2c<-glm(H.wet~log.lake.area, family=gaussian, data=wet.mod.df)
m3c<-glm(H.wet~max.depth, family=gaussian, data=wet.mod.df)
m4c<-glm(H.wet~IWS.stream.density, family=gaussian, data=wet.mod.df)
m5c<-glm(H.wet~log.area.upstream.lakes, family=gaussian, data=wet.mod.df)
m6c<-glm(H.wet~perc.agric.dev.land, family=gaussian, data=wet.mod.df)
m7c<-glm(H.wet~julian.day, family=gaussian, data=wet.mod.df)
m8c<-glm(H.wet~ log.lake.area+max.depth+IWS.stream.density+log.area.upstream.lakes+perc.agric.dev.land+julian.day+year, family = gaussian, data = wet.mod.df)
m9c<-glm(H.wet~1,family=gaussian, data=wet.mod.df)

# Model selection
output.wet.H<-model.sel(m1c,m2c,m3c,m4c,m5c,m6c,m7c,m8c,m9c)

# Likelihood ratio test
lrtest(m7c,m9c)

# Candidate models: Species richness - Wetland species
m1d<-glm.nb(richness.wet~year, data=wet.mod.df)
m2d<-glm.nb(richness.wet~log.lake.area, data=wet.mod.df)
m3d<-glm.nb(richness.wet~max.depth, data=wet.mod.df)
m4d<-glm.nb(richness.wet~IWS.stream.density, data=wet.mod.df)
m5d<-glm.nb(richness.wet~log.area.upstream.lakes, data=wet.mod.df)
m6d<-glm.nb(richness.wet~perc.agric.dev.land, data=wet.mod.df)
m7d<-glm.nb(richness.wet~julian.day, data=wet.mod.df)
m8d<-glm.nb(richness.wet~ log.lake.area+max.depth+IWS.stream.density+log.area.upstream.lakes+perc.agric.dev.land+julian.day+year, data = wet.mod.df)
m9d<-glm.nb(richness.wet~1, data=wet.mod.df)

# Model selection
output.wet.rich<-model.sel(m1d,m2d,m3d,m4d,m5d,m6d,m7d,m8d,m9d)

# Likelihood ratio test
lrtest(m7d,m9d)


##################################
#     3.2 AQUATIC MODELS         #
##################################

# Transform predictors
aqu.mod.df$log.lake.area<-log(aqu.mod.df$lake.area)
aqu.mod.df$log.area.upstream.lakes<-log(aqu.mod.df$area.upstream.lakes+1)

# Ensure correct data types
aqu.mod.df$H.aqu<-as.numeric(aqu.mod.df$H.aqu)
aqu.mod.df$richness.aqu<-as.integer(aqu.mod.df$richness.aqu)
aqu.mod.df$year<-as.factor(aqu.mod.df$year)

# Full models (used to assess multicollinearity)
full.aqu.H<-glm(H.aqu ~ log.lake.area + max.depth + IWS.stream.density + log.area.upstream.lakes + 
              perc.agric.dev.land + julian.day, family = gaussian, data = aqu.mod.df)

vif(full.aqu.H)

full.aqu.rich<-glm.nb(richness.aqu~ log.lake.area+max.depth+IWS.stream.density+log.area.upstream.lakes+
                 perc.agric.dev.land+julian.day, data = aqu.mod.df)
vif(full.aqu.rich)


# Candidate models: Shannon Diversity - Aquatic species

m1e<-glm(H.aqu~year, family=gaussian, data=aqu.mod.df)
m2e<-glm(H.aqu~log.lake.area, family=gaussian, data=aqu.mod.df)
m3e<-glm(H.aqu~max.depth, family=gaussian, data=aqu.mod.df)
m4e<-glm(H.aqu~IWS.stream.density, family=gaussian, data=aqu.mod.df)
m5e<-glm(H.aqu~log.area.upstream.lakes, family=gaussian, data=aqu.mod.df)
m6e<-glm(H.aqu~perc.agric.dev.land, family=gaussian, data=aqu.mod.df)
m7e<-glm(H.aqu~julian.day, family=gaussian, data=aqu.mod.df)
m8e<-glm(H.aqu~ log.lake.area+max.depth+IWS.stream.density+log.area.upstream.lakes+perc.agric.dev.land+julian.day+year, family = gaussian, data = aqu.mod.df)
m9e<-glm(H.aqu~1,family=gaussian, data=aqu.mod.df)

# Model selection
output.aqu.H<-model.sel(m1e,m2e,m3e,m4e,m5e,m6e,m7e,m8e,m9e)

# Likelihood ratio test
lrtest(m7e,m9e)
lrtest(m2e,m9e)
lrtest(m3e,m9e)
lrtest(m5e,m9e)


# Candidate models: Species richness - Aquatic species
m1f<-glm.nb(richness.aqu~year, data=aqu.mod.df)
m2f<-glm.nb(richness.aqu~log.lake.area, data=aqu.mod.df)
m3f<-glm.nb(richness.aqu~max.depth, data=aqu.mod.df)
m4f<-glm.nb(richness.aqu~IWS.stream.density, data=aqu.mod.df)
m5f<-glm.nb(richness.aqu~log.area.upstream.lakes, data=aqu.mod.df)
m6f<-glm.nb(richness.aqu~perc.agric.dev.land, data=aqu.mod.df)
m7f<-glm.nb(richness.aqu~julian.day, data=aqu.mod.df)
m8f<-glm.nb(richness.aqu~ log.lake.area+max.depth+IWS.stream.density+log.area.upstream.lakes+perc.agric.dev.land+julian.day+year, data = aqu.mod.df)
m9f<-glm.nb(richness.aqu~1, data=aqu.mod.df)

# Model selection
output.aqu.rich<-model.sel(m1f,m2f,m3f,m4f,m5f,m6f,m7f,m8f,m9f)

# Likelihood ratio test
lrtest(m2f,m9f)
lrtest(m7f,m9f)
lrtest(m3f,m9f)
lrtest(m5f,m9f)

#################################
# Plots Models Alpha Diversity  #
#################################


############################################
# COMMON THEME (apply to all plots)       #
############################################

theme_alpha <- theme(
  legend.title = element_blank(),
  panel.background = element_rect(fill = "white"),
  panel.border = element_rect(linetype = "solid", fill = NA),
  axis.text = element_text(size = 16),
  axis.title.y = element_text(margin = margin(r = 20), size = 24),
  axis.title.x = element_text(margin = margin(t = 20), size = 24),
  plot.margin = margin(t = 20, r = 20, b = 20, l = 20))

##################################
#  Plot Terrestrial Communities  #
##################################

###########################
#    Shannon-Diversity    #
###########################

m7a.pred <- ggpredict(m7a, terms = c("julian.day [211,222,198,154,180,140,234,148,242,161,226,137,144,221,227,144,169,137,214,219,175,218]"))

julian.day.terr.shannon <- ggplot(m7a.pred, aes(x = x, y = predicted)) +
  geom_line(linewidth = 1.5, color = "darkgray") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  stat_summary(fun = mean, fill = "coral", color = "black", geom = "point", size = 5, shape = 21) +
  geom_point(data = terr.mod.df,
             aes(x = julian.day, y = H.terr),
             size = 5, color = "black", shape = 21, stroke = 0.8) +
  ylab("Shannon diversity") +
  xlab("Julian day") +
  ylim(0, 3) +
  theme_alpha

#############################
#     Species Richness      #
#############################

m7b.pred <- ggpredict(m7b, terms = c("julian.day [211,222,198,154,180,140,234,148,242,161,226,137,144,221,227,144,169,137,214,219,175,218]"))

julian.day.terr.richness <- ggplot(m7b.pred, aes(x = x, y = predicted)) +
  geom_line(linewidth = 1.5, color = "darkgray") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  stat_summary(fun = mean, fill = "coral", color = "black", geom = "point", size = 5, shape = 21) +
  geom_point(data = terr.mod.df,
             aes(x = julian.day, y = richness.terr),
             size = 5, color = "black", shape = 21, stroke = 0.8) +
  ylab("Species richness") +
  xlab("Julian day") +
  ylim(0, 60) +
  theme_alpha +
  theme(axis.title.x = element_blank())

##################################
#    Plot Wetland Communities    #
##################################

#############################
#     Shannon-Diversity     #
#############################

m7c.pred <- ggpredict(m7c, terms = c("julian.day [211,222,198,154,180,140,234,148,242,161,226,137,144,221,227,144,169,137,214,219,175,218]"))

julian.day.wet.shannon <- ggplot(m7c.pred, aes(x = x, y = predicted)) +
  geom_line(linewidth = 1.5, color = "darkgray") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  stat_summary(fun = mean, fill = "coral", color = "black", geom = "point", size = 5, shape = 21) +
  geom_point(data = wet.mod.df,
             aes(x = julian.day, y = H.wet),
             size = 5, color = "black", shape = 21, stroke = 0.8) +
  ylab("Shannon diversity") +
  xlab("Julian day") +
  ylim(0, 3) +
  theme_alpha +
  theme(axis.title.y = element_blank())


#############################
#      Species Richness     #
#############################

m7d.pred <- ggpredict(m7d, terms = c("julian.day [211,222,198,154,180,140,234,148,242,161,226,137,144,221,227,144,169,137,214,219,175,218]"))

julian.day.wet.richness <- ggplot(m7d.pred, aes(x = x, y = predicted)) +
  geom_line(linewidth = 1.5, color = "darkgray") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  stat_summary(fun = mean, fill = "coral", color = "black", geom = "point", size = 5, shape = 21) +
  geom_point(data = wet.mod.df,
             aes(x = julian.day, y = richness.wet),
             size = 5, color = "black", shape = 21, stroke = 0.8) +
  ylab("Species richness") +
  xlab("Julian day") +
  ylim(0, 60) +
  theme_alpha +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())


#################################
# Combine plots                #
#################################

combinedPlot <- plot_grid(
  julian.day.terr.richness,
  julian.day.wet.richness,
  julian.day.terr.shannon,
  julian.day.wet.shannon,  
  labels = "auto",
  label_size = 25,
  ncol = 2,
  align = "hv"
)

png("Figure 3 Alpha Diversity Models.png",
    width = 4000, height = 3000, units = "px", res = 300)

print(combinedPlot)

dev.off()

