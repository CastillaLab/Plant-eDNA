# First, load all packages
library(tidyverse)
library(dplyr)
library(vegan)
library(MuMIn)
library(lmerTest)
library(lmtest)
library(ggeffects)
library(ggplot2)
library(cowplot)
library(car)
library(MASS)

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

#####################################################
# SUBSET TERRESTRIAL, AQUATIC, AND WETLAND SPECIES  #
#####################################################

# Load the list with the species, their life habit, and status
spinfo<-read.csv("spinfo.csv", header=TRUE)

#######################
# TERRESTRIAL SPECIES #
#######################

# We create a list with the terrestrial species 
keep_spp_terr<-spinfo$species[spinfo$habitat == "terrestrial"]
length(keep_spp_terr)

# We keep only the terrestrial species in our plant community matrix
CM.filt.terr<-CM.filtered[,c(1:3,which(colnames(CM.filtered) %in% keep_spp_terr))]
head(CM.filt.terr)

#######################
#   WETLAND SPECIES   #
#######################

# We create a list with the wetland species 
keep_spp_wet<-spinfo$species[spinfo$habitat == "wetland"]
length(keep_spp_wet)

# We keep only the wetland species in our plant community matrix
CM.filt.wet<-CM.filtered[,c(1:3,which(colnames(CM.filtered) %in% keep_spp_wet))]
head(CM.filt.wet)


#######################
#   AQUATIC SPECIES   #
#######################

# We create a list with the aquatic species 
keep_spp_aqu<-spinfo$species[spinfo$habitat == "aquatic"]
length(keep_spp_aqu)

# We keep only the wetland species in our plant community matrix
CM.filt.aqu<-CM.filtered[,c(1:3,which(colnames(CM.filtered) %in% keep_spp_aqu))]
head(CM.filt.aqu)


######################################
##     Summarize for each lake       #
######################################

# Terrestrial plants
CM.filt.terr<-CM.filt.terr[,c(1,4:134)]
CM.filt.terr.lake<-CM.filt.terr %>% group_by(lake_name) %>% summarise(across(everything(),list(sum)))


# Wetland plants
CM.filt.wet<-CM.filt.wet[,c(1,4:63)]
CM.filt.wet.lake<-CM.filt.wet %>% group_by(lake_name) %>% summarise(across(everything(),list(sum)))


# Aquatic plants
CM.filt.aqu<-CM.filt.aqu[,c(1,4:33)]
CM.filt.aqu.lake<-CM.filt.aqu %>% group_by(lake_name) %>% summarise(across(everything(),list(sum)))


#####################################
##   2. CALCULATE ALPHA DIVERSITY  ##  
#####################################

################################
## 2.1 Terrestrial community  ##
################################

#CM.filt.terr.lake<-read.csv("CM.filt.terr.lake.csv", header=TRUE)

df.terr<-CM.filt.terr.lake[,2:132]

# Shannon's H'
H.terr <- diversity(df.terr)

# Observed Richness
richness.terr <- specnumber(df.terr)  

lake_name<-CM.filt.terr.lake$lake_name

merg.diversity.terr<-cbind(lake_name,H.terr,richness.terr)

lake_traits<-read.csv("lakes traits.csv", header = T)

terr.mod.df<-merge(lake_traits,merg.diversity.terr,by="lake_name")


################################
##   2.2 Wetland community    ##
################################

#CM.filt.wet.lake<-read.csv("CM.filt.wet.lake.csv", header=TRUE)

df.wet<-CM.filt.wet.lake[,2:61]

# Shannon's H'
H.wet <- diversity(df.wet)

# Observed Richness
richness.wet <- specnumber(df.wet)  

merg.diversity.wet<-cbind(lake_name,H.wet,richness.wet)

wet.mod.df<-merge(lake_traits,merg.diversity.wet,by="lake_name")

################################
##   2.3 Aquatic community    ##
################################

#CM.filt.aqu.lake<-read.csv("CM.filt.aqu.lake.csv", header=TRUE)

df.aqu<-CM.filt.aqu.lake[,2:31]

# Shannon's H'
H.aqu <- diversity(df.aqu)

# Observed Richness
richness.aqu <- specnumber(df.aqu)  

merg.diversity.aqu<-cbind(lake_name,H.aqu,richness.aqu)

aqu.mod.df<-merge(lake_traits,merg.diversity.aqu,by="lake_name")


#####################################
##  3. MODELS ALPHA DIVERSITY      ##  
#####################################

##################################
#         3.1 Terrestrial        #
##################################

terr.mod.df$log.lake.area<-log(terr.mod.df$lake.area)
terr.mod.df$log.area.upstream.lakes<-log(terr.mod.df$area.upstream.lakes+1)

terr.mod.df$H.terr<-as.numeric(terr.mod.df$H.terr)
terr.mod.df$richness.terr<-as.integer(terr.mod.df$richness.terr)

full.terr.H<-glm(H.terr ~ log.lake.area + max.depth + IWS.stream.density + log.area.upstream.lakes + 
              perc.agric.dev.land + julian.day, family = gaussian, data = terr.mod.df)

vif(full.terr.H)

full.terr.rich<-glm.nb(richness.terr~ log.lake.area+max.depth+IWS.stream.density+log.area.upstream.lakes+
                 perc.agric.dev.land+julian.day, data = terr.mod.df)
vif(full.terr.rich)


#Shannon-Wiener diversity index
m1a<-glm(H.terr~log.lake.area, family=gaussian, data=terr.mod.df)
m2a<-glm(H.terr~max.depth, family=gaussian, data=terr.mod.df)
m3a<-glm(H.terr~IWS.stream.density, family=gaussian, data=terr.mod.df)
m4a<-glm(H.terr~log.area.upstream.lakes, family=gaussian, data=terr.mod.df)
m5a<-glm(H.terr~perc.agric.dev.land, family=gaussian, data=terr.mod.df)
m6a<-glm(H.terr~julian.day, family=gaussian, data=terr.mod.df)
m7a<-glm(H.terr~ log.lake.area+max.depth+IWS.stream.density+log.area.upstream.lakes+perc.agric.dev.land+julian.day, family = gaussian, data = terr.mod.df)
m8a<-glm(H.terr~1,family=gaussian, data=terr.mod.df)

output.terr.H<-model.sel(m1a,m2a,m3a,m4a,m5a,m6a,m7a,m8a)

output.terr.H

lrtest(m6a,m8a)


#############################
# Plot Predictions Shannon  #
#############################

# Extract the predicted values
m6a.pred<-ggpredict(m6a, terms = c("julian.day [211,222,198,154,180,140,234,148,242,161,226,137,144,221,227,144,169,137,214,219,175,218]"))

# Save the predicted values and add the observed values 
write.csv(m6a.pred,"julian.day.pred.terr.csv",row.names=FALSE)

# Load the data with predicted and observed values for the response variable
m6a.pred.obs<-read.csv("julian.day.pred.terr.csv", header=TRUE)


#julian.day

julian.day.shannon.terr<-ggplot(m6a.pred.obs, aes(x=x, y=predicted))+
  geom_line(linewidth=1.5, color="darkgray")+geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=0.2)+
  stat_summary(fun=mean, fill="coral", color= "black", geom="point", size=5, shape=21)+
  geom_point(aes(x=x,y=H.terr),size=5,color="black",shape=21,stroke=0.8)+
  theme(legend.title = element_blank(),panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", fill = NA))+
  ylab("Shannon diversity")+
  xlab("Julian day")+
  theme(axis.title = element_text(face="bold", size=24),
        axis.text = element_text(size=16))

julian.day.shannon.terr

ggsave("Shannon_terrestrial community.png", width=1600, height=1200, units = "px", dpi = 300)



#Species richness

m1b<-glm.nb(richness.terr~log.lake.area, data=terr.mod.df)
m2b<-glm.nb(richness.terr~max.depth, data=terr.mod.df)
m3b<-glm.nb(richness.terr~IWS.stream.density, data=terr.mod.df)
m4b<-glm.nb(richness.terr~log.area.upstream.lakes, data=terr.mod.df)
m5b<-glm.nb(richness.terr~perc.agric.dev.land, data=terr.mod.df)
m6b<-glm.nb(richness.terr~julian.day, data=terr.mod.df)
m7b<-glm.nb(richness.terr~ log.lake.area+max.depth+IWS.stream.density+log.area.upstream.lakes+perc.agric.dev.land+julian.day, data = terr.mod.df)
m8b<-glm.nb(richness.terr~1, data=terr.mod.df)

output.terr.rich<-model.sel(m1b,m2b,m3b,m4b,m5b,m6b,m7b,m8b)

output.terr.rich

lrtest(m6b,m8b)

#############################
# Plot Predictions Richness #
#############################

# Extract the predicted values
m6b.pred<-ggpredict(m6b, terms = c("julian.day [211,222,198,154,180,140,234,148,242,161,226,137,144,221,227,144,169,137,214,219,175,218]"))

# Save the predicted values and add the observed values
write.csv(m6b.pred,"julian.day.pred.richness.terr.csv",row.names=FALSE)

# Load the data with predicted and observed values for the response variable
m6b.pred.obs<-read.csv("julian.day.pred.richness.terr.csv", header=TRUE)



julian.day.terr.richness<-ggplot(m6d.pred.obs, aes(x=x, y=predicted))+
  geom_line(size=1.5, color="darkgray")+geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=0.2)+
  stat_summary(fun=mean, fill="coral", color= "black", geom="point", size=5, shape=21)+
  geom_point(aes(x=x,y=richness.terr),size=5,color="black",shape=21,stroke=0.8)+
  theme(legend.title = element_blank(),panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", fill = NA))+
  ylab("Species richness")+
  xlab("Julian day")+
  theme(axis.title = element_text(face="bold", size=24),
        axis.text = element_text(size=16))

julian.day.terr.richness

ggsave("Richness_terrestrial community.png", width=1600, height=1200, units = "px", dpi = 300)


##################################
#           3.2 Wetland          #
##################################

wet.mod.df$log.lake.area<-log(wet.mod.df$lake.area)
wet.mod.df$log.area.upstream.lakes<-log(wet.mod.df$area.upstream.lakes+1)

wet.mod.df$H.wet<-as.numeric(wet.mod.df$H.wet)
wet.mod.df$richness.wet<-as.integer(wet.mod.df$richness.wet)

full.wet.H<-glm(H.wet ~ log.lake.area + max.depth + IWS.stream.density + log.area.upstream.lakes + 
              perc.agric.dev.land + julian.day, family = gaussian, data = wet.mod.df)

vif(full.wet.H)

full.wet.rich<-glm.nb(richness.wet~ log.lake.area+max.depth+IWS.stream.density+log.area.upstream.lakes+
                 perc.agric.dev.land+julian.day, data = wet.mod.df)
vif(full.wet.rich)



# Shannon diversity index

m1c<-glm(H.wet~log.lake.area, family=gaussian, data=wet.mod.df)
m2c<-glm(H.wet~max.depth, family=gaussian, data=wet.mod.df)
m3c<-glm(H.wet~IWS.stream.density, family=gaussian, data=wet.mod.df)
m4c<-glm(H.wet~log.area.upstream.lakes, family=gaussian, data=wet.mod.df)
m5c<-glm(H.wet~perc.agric.dev.land, family=gaussian, data=wet.mod.df)
m6c<-glm(H.wet~julian.day, family=gaussian, data=wet.mod.df)
m7c<-glm(H.wet~ log.lake.area+max.depth+IWS.stream.density+log.area.upstream.lakes+perc.agric.dev.land+julian.day, family = gaussian, data = wet.mod.df)
m8c<-glm(H.wet~1,family=gaussian, data=wet.mod.df)


output.wet.H<-model.sel(m1c,m2c,m3c,m4c,m5c,m6c,m7c,m8c)
output.wet.H

lrtest(m6c,m8c)

#############################
# Plot Predictions Shannon  #
#############################

# Extract the predicted values
m6c.pred<-ggpredict(m6c, terms = c("julian.day [211,222,198,154,180,140,234,148,242,161,226,137,144,221,227,144,169,137,214,219,175,218]"))

# Save the predicted values and add the observed values
write.csv(m6c.pred,"julian.day.pred.wet.csv",row.names=FALSE)

# Load the data with predicted and observed values for the response variable
m6c.pred.obs<-read.csv("julian.day.pred.wet.csv", header=TRUE)


#julian.day

julian.day.shannon.wet<-ggplot(m6c.pred.obs, aes(x=x, y=predicted))+
  geom_line(linewidth=1.5, color="darkgray")+geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=0.2)+
  stat_summary(fun=mean, fill="coral", color= "black", geom="point", size=5, shape=21)+
  geom_point(aes(x=x,y=H.wet),size=5,color="black",shape=21,stroke=0.8)+
  theme(legend.title = element_blank(),panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", fill = NA))+
  ylab("Shannon diversity")+
  xlab("Julian day")+
  theme(axis.title = element_text(face="bold", size=24),
        axis.text = element_text(size=16))

julian.day.shannon.wet

ggsave("Shannon_wetland community.png", width=1600, height=1200, units = "px", dpi = 300)


# Species richness

m1d<-glm.nb(richness.wet~log.lake.area, data=wet.mod.df)
m2d<-glm.nb(richness.wet~max.depth, data=wet.mod.df)
m3d<-glm.nb(richness.wet~IWS.stream.density, data=wet.mod.df)
m4d<-glm.nb(richness.wet~log.area.upstream.lakes, data=wet.mod.df)
m5d<-glm.nb(richness.wet~perc.agric.dev.land, data=wet.mod.df)
m6d<-glm.nb(richness.wet~julian.day, data=wet.mod.df)
m7d<-glm.nb(richness.wet~ log.lake.area+max.depth+IWS.stream.density+log.area.upstream.lakes+perc.agric.dev.land+julian.day, data = wet.mod.df)
m8d<-glm.nb(richness.wet~1, data=wet.mod.df)

output.wet.rich<-model.sel(m1d,m2d,m3d,m4d,m5d,m6d,m7d,m8d)

output.wet.rich

lrtest(m6d,m8d)

#############################
# Plot Predictions Richness #
#############################

# Extract the predicted values
m6d.pred<-ggpredict(m6d, terms = c("julian.day [211,222,198,154,180,140,234,148,242,161,226,137,144,221,227,144,169,137,214,219,175,218]"))

# Save the predicted values and add the observed values
write.csv(m6d.pred,"julian.day.pred.richness.wet.csv",row.names=FALSE)

# Load the data with predicted and observed values for the response variable
m6d.pred.obs<-read.csv("julian.day.pred.richness.wet.csv", header=TRUE)



julian.day.wet.richness<-ggplot(m6d.pred.obs, aes(x=x, y=predicted))+
  geom_line(size=1.5, color="darkgray")+geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=0.2)+
  stat_summary(fun=mean, fill="coral", color= "black", geom="point", size=5, shape=21)+
  geom_point(aes(x=x,y=richness.wet),size=5,color="black",shape=21,stroke=0.8)+
  theme(legend.title = element_blank(),panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", fill = NA))+
  ylab("Species richness")+
  xlab("Julian day")+
  theme(axis.title = element_text(face="bold", size=24),
        axis.text = element_text(size=16))

julian.day.wet.richness

ggsave("Richness_wetland community.png", width=1600, height=1200, units = "px", dpi = 300)



##################################
#           3.3 Aquatic          #
##################################

aqu.mod.df$log.lake.area<-log(aqu.mod.df$lake.area)
aqu.mod.df$log.area.upstream.lakes<-log(aqu.mod.df$area.upstream.lakes+1)

aqu.mod.df$H.aqu<-as.numeric(aqu.mod.df$H.aqu)
aqu.mod.df$richness.aqu<-as.integer(aqu.mod.df$richness.aqu)

full.aqu.H<-glm(H.aqu ~ log.lake.area + max.depth + IWS.stream.density + log.area.upstream.lakes + 
              perc.agric.dev.land + julian.day, family = gaussian, data = aqu.mod.df)

vif(full.aqu.H)

full.aqu.rich<-glm.nb(richness.aqu~ log.lake.area+max.depth+IWS.stream.density+log.area.upstream.lakes+
                 perc.agric.dev.land+julian.day, data = aqu.mod.df)
vif(full.aqu.rich)


# Shannon diversity index

m1e<-glm(H.aqu~log.lake.area, family=gaussian, data=aqu.mod.df)
m2e<-glm(H.aqu~max.depth, family=gaussian, data=aqu.mod.df)
m3e<-glm(H.aqu~IWS.stream.density, family=gaussian, data=aqu.mod.df)
m4e<-glm(H.aqu~log.area.upstream.lakes, family=gaussian, data=aqu.mod.df)
m5e<-glm(H.aqu~perc.agric.dev.land, family=gaussian, data=aqu.mod.df)
m6e<-glm(H.aqu~julian.day, family=gaussian, data=aqu.mod.df)
m7e<-glm(H.aqu~ log.lake.area+max.depth+IWS.stream.density+log.area.upstream.lakes+perc.agric.dev.land+julian.day, family = gaussian, data = aqu.mod.df)
m8e<-glm(H.aqu~1,family=gaussian, data=aqu.mod.df)

output.aqu.H<-model.sel(m1e,m2e,m3e,m4e,m5e,m6e,m7e,m8e)
output.aqu.H

lrtest(m1e,m8e)
lrtest(m2e,m8e)
lrtest(m4e,m8e)
lrtest(m6e,m8e)

# Species richness

m1f<-glm.nb(richness.aqu~log.lake.area, data=aqu.mod.df)
m2f<-glm.nb(richness.aqu~max.depth, data=aqu.mod.df)
m3f<-glm.nb(richness.aqu~IWS.stream.density, data=aqu.mod.df)
m4f<-glm.nb(richness.aqu~log.area.upstream.lakes, data=aqu.mod.df)
m5f<-glm.nb(richness.aqu~perc.agric.dev.land, data=aqu.mod.df)
m6f<-glm.nb(richness.aqu~julian.day, data=aqu.mod.df)
m7f<-glm.nb(richness.aqu~ log.lake.area+max.depth+IWS.stream.density+log.area.upstream.lakes+perc.agric.dev.land+julian.day, data = aqu.mod.df)
m8f<-glm.nb(richness.aqu~1, data=aqu.mod.df)

output.aqu.rich<-model.sel(m1f,m2f,m3f,m4f,m5f,m6f,m7f,m8f)

output.aqu.rich

lrtest(m1f,m8f)
lrtest(m6f,m8f)


##############################
# Plot Predictions Richness  #
##############################

# Extract the predicted values
m6f.pred<-ggpredict(m6f, terms = c("julian.day [211,222,198,154,180,140,234,148,242,161,226,137,144,221,227,144,169,137,214,219,175,218]"))
m1f.pred<-ggpredict(m1f, terms = c("log.lake.area [6.100274,7.464928,6.274235,4.580775,4.500254,2.953347,3.958907,2.700690,8.325214,6.186476,9.002903,4.071417,4.456902,8.339530,8.819115,3.925926,5.288671,4.821329,4.671426,8.933135,7.524788,5.508376]"))

# Save the predicted values and add the observed values
write.csv(m6f.pred,"julian.day.pred.aquatic.csv",row.names=FALSE)
write.csv(m1f.pred,"lake.area.pred.aquatic.csv",row.names=FALSE)

# Load the data with predicted and observed values for the response variable
m6f.pred.obs<-read.csv("julian.day.pred.aquatic.csv", header=TRUE)
m1f.pred.obs<-read.csv("lake.area.pred.aquatic.csv", header=TRUE)


#julian.day

julian.day.richness.aqu<-ggplot(m6f.pred.obs, aes(x=x, y=predicted))+
  geom_line(linewidth=1.5, color="darkgray")+geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=0.2)+
  stat_summary(fun=mean, fill="coral", color= "black", geom="point", size=5, shape=21)+
  geom_point(aes(x=x,y=richness.aqu),size=5,color="black",shape=21,stroke=0.8)+
  theme(legend.title = element_blank(),panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", fill = NA))+
  ylab("Spcies richness")+
  xlab("Julian day")+
  theme(axis.title = element_text(face="bold", size=24),
        axis.text = element_text(size=16))

julian.day.richness.aqu


#log.lake.area

lake.area.richness.aqu<-ggplot(m1f.pred.obs, aes(x=x, y=predicted))+
  geom_line(size=1.5, color="darkgray")+geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=0.2)+
  stat_summary(fun=mean, fill="coral", color= "black", geom="point", size=5, shape=21)+
  geom_point(aes(x=x,y=richness.aqu),size=5,color="black",shape=21,stroke=0.8)+
  theme(legend.title = element_blank(),panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", fill = NA))+
  xlab("log(lake area)")+
  theme(axis.title.x = element_text(face="bold", size=24),
        axis.title.y = element_blank(),
        axis.text = element_text(size=16))

lake.area.richness.aqu

combinedPlot2 <- plot_grid(julian.day.richness.aqu,lake.area.richness.aqu, labels='auto', label_y= 1, label_x=0.11, align='hv', label_size=20, 
                           ncol=2, rel_widths=1, rel_heights= 1)

png("Richnes_aquatic community.png",width = 4000, height = 1200, units = "px", res=300)
combinedPlot2
dev.off()
