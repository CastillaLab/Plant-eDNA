#Load packages
library(dplyr)
library(tidyverse)
library(ggplot2)
library(GGally)
library(MuMIn)
library(lmerTest)
library(lmtest)
library(AER)
library(MASS)
library(car)
library(ggeffects)
library(ggplot2)
library(cowplot)
library(car)


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

# Terrestrial plants
CM.filt.terr<-CM.filt.terr[,c(1,4:131)]
CM.filt.terr.lake<-CM.filt.terr %>% group_by(lake_name) %>% summarise(across(everything(),list(sum)))

# Wetland plants
CM.filt.wet<-CM.filt.wet[,c(1,4:61)]
CM.filt.wet.lake<-CM.filt.wet %>% group_by(lake_name) %>% summarise(across(everything(),list(sum)))

# Aquatic plants
CM.filt.aqu<-CM.filt.aqu[,c(1,4:33)]
CM.filt.aqu.lake<-CM.filt.aqu %>% group_by(lake_name) %>% summarise(across(everything(),list(sum)))


# For each lake, we sum the sequences across all columns so we obtain the total 
# number of sequences per lake. We do this for each dataset, so we obtained the
# total number of plant species, terrestrial species, wetland species and
# aquatic species for each lake.

CM.filt.terr.lake$terrestrial.reads=rowSums(CM.filt.terr.lake[,2:129])
CM.filt.wet.lake$wetland.reads=rowSums(CM.filt.wet.lake[,2:59])
CM.filt.aqu.lake$aquatic.reads=rowSums(CM.filt.aqu.lake[,2:31])

# Create a data frame with the total number of reads per lake, the total number 
# of reads for terrestrial plants, the total number of reads for wetland plants,
# the total number of reads for aquatic plants, the proportion of reads of 
# terrestrial plants, the proportion of reads of wetland plants and the 
# proportion of reads of aquatic plants.

df_list<-list(CM.filt.terr.lake,CM.filt.wet.lake,CM.filt.aqu.lake)
df.prop<-df_list %>% reduce(full_join, by='lake_name')
df.prop<-df.prop %>% select(lake_name,terrestrial.reads,wetland.reads,aquatic.reads)

df.prop$total.reads<-df.prop$terrestrial.reads+df.prop$wetland.reads+df.prop$aquatic.reads

df.prop$prop.terrestrial<-df.prop$terrestrial.reads/df.prop$total.reads
df.prop$prop.wetland<-df.prop$wetland.reads/df.prop$total.reads
df.prop$prop.aquatic<-df.prop$aquatic.reads/df.prop$total.reads

# Add the lake traits and landscape traits to the dataframe with the 
# proportions

lake.traits<-read.csv("lakes traits.csv", header=TRUE)
merged.df<-merge(df.prop,lake.traits,by="lake_name")

############################################
# Model selection for terrestrial species  #
############################################

merged.df$forest.500<-merged.df$perc_forest_2011_buffer500
merged.df$developed.500<-merged.df$perc_dev_land_2011_buffer500
merged.df$agric.500<-merged.df$perc_agric_land_2011_buffer500

merged.df$log.lake.area<-log(merged.df$lake.area)


# VIF 

full.terr<-glm(prop.terrestrial~log.lake.area+max.depth+forest.500+
            agric.500+developed.500,family=gaussian, data=merged.df)

vif(full.terr)

## Forest percentage

m1a<-glm(prop.terrestrial~log.lake.area, family=gaussian, data=merged.df)
m2a<-glm(prop.terrestrial~max.depth, family=gaussian, data=merged.df)
m3a<-glm(prop.terrestrial~forest.500, family=gaussian, data=merged.df)
m4a<-glm(prop.terrestrial~developed.500, family=gaussian, data=merged.df)
m5a<-glm(prop.terrestrial~agric.500, family=gaussian, data=merged.df)
m6a<-glm(prop.terrestrial~log.lake.area+max.depth+forest.500, family=gaussian, data = merged.df)
m7a<-glm(prop.terrestrial~agric.500+developed.500, family=gaussian, data = merged.df) 
m8a<-glm(prop.terrestrial~log.lake.area+max.depth+forest.500+agric.500+developed.500, family=gaussian, data=merged.df)
m9a<-glm(prop.terrestrial~1, family=gaussian, data=merged.df)

output.terr<-model.sel(m1a,m2a,m3a,m4a,m5a,m6a,m7a,m8a,m9a)

output.terr

lrtest(m4a,m9a)
lrtest(m3a,m9a)


############################################
#   Model selection for wetland species    #
############################################


# VIF 

full.wet<-glm(prop.wetland~log.lake.area+max.depth+forest.500+
            agric.500+developed.500,family=gaussian, data=merged.df)

vif(full.wet)


m1b<-glm(prop.wetland~log.lake.area, family=gaussian, data=merged.df)
m2b<-glm(prop.wetland~max.depth, family=gaussian, data=merged.df)
m3b<-glm(prop.wetland~forest.500, family=gaussian, data=merged.df)
m4b<-glm(prop.wetland~developed.500, family=gaussian, data=merged.df)
m5b<-glm(prop.wetland~agric.500, family=gaussian, data=merged.df)
m6b<-glm(prop.wetland~log.lake.area+max.depth+forest.500,
         family=gaussian, data = merged.df)
m7b<-glm(prop.wetland~agric.500+developed.500,
         family=gaussian, data = merged.df) 
m8b<-glm(prop.wetland~log.lake.area+max.depth+forest.500+
           agric.500+developed.500, family=gaussian, data=merged.df)
m9b<-glm(prop.wetland~1, family=gaussian, data=merged.df)

output.wet<-model.sel(m1b,m2b,m3b,m4b,m5b,m6b,m7b,m8b,m9b)

output.wet

lrtest(m3b,m9b)


############################################
#   Model selection for aquatic species    #
############################################

# VIF 

full.aqu<-glm(prop.aquatic~log.lake.area+max.depth+forest.500+
            agric.500+developed.500,family=gaussian, data=merged.df)

vif(full.aqu)

## Forest percentage

m1c<-glm(prop.aquatic~log.lake.area, family=gaussian, data=merged.df)
m2c<-glm(prop.aquatic~max.depth, family=gaussian, data=merged.df)
m3c<-glm(prop.aquatic~forest.500, family=gaussian, data=merged.df)
m4c<-glm(prop.aquatic~developed.500, family=gaussian, data=merged.df)
m5c<-glm(prop.aquatic~agric.500, family=gaussian, data=merged.df)
m6c<-glm(prop.aquatic~log.lake.area+max.depth+forest.500, family=gaussian,
         data = merged.df)
m7c<-glm(prop.aquatic~agric.500+developed.500, family=gaussian,
         data = merged.df) 
m8c<-glm(prop.aquatic~log.lake.area+max.depth+forest.500+agric.500+
           developed.500,family=gaussian, data=merged.df)
m9c<-glm(prop.aquatic~1, family=gaussian, data=merged.df)

output.aqu<-model.sel(m1c,m2c,m3c,m4c,m5c,m6c,m7c,m8c,m9c)

output.aqu

lrtest(m4c,m9c)
lrtest(m7c,m4c)


#############
#   PLOTS   #
#############

#########################
#  Terrestrial species  #
#########################

# Developed land

m4a.pred<-ggpredict(m4a, terms = c("developed.500 [46.26,6.33,72.81,7.39,9.08,17.68,24.04,81.25,33.09,22.95,49.12,7.2,20.18,10.68,19.82,9.15,28.32,13.98,92.2,21.65,12.42,1.5]]"))
write.csv(m4a.pred,"terr.dev.all.pred.csv",row.names=FALSE)
m4a.pred.obs<-read.csv("terr.dev.all.pred.csv", header=TRUE)

p1<-ggplot(m4a.pred.obs, aes(x=x, y=predicted))+
  geom_line(linewidth=1.5, color="darkgray")+geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=0.2)+
  stat_summary(fun=mean, fill="coral", color= "black", geom="point", size=5, shape=21)+
  geom_point(aes(x=x,y=prop.terrestrial), size=5,shape=21,stroke=0.8)+
  ylim(0,1)+
  theme(legend.title = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.title = element_text(size=16),
        axis.text = element_text(size=14),
        panel.border = element_rect(linetype = "solid", fill = NA))+
  ylab("Proportion of \n terrestrial species")+
  xlab("% Developed land (500 m)")

# Forest cover

m3a.pred<-ggpredict(m3a, terms = c("forest.500 [18.76,44,6.38,20.33,46.93,3.58,57.25,5.55,52.87,39.92,9.99,49.31,26.87,34.1,37.03,57.94,37.6,58.82,1.19,52.4,58.69,28.55]"))
write.csv(m3a.pred,"terr.forest.all.pred.csv",row.names=FALSE)
m3a.pred.obs<-read.csv("terr.forest.all.pred.csv", header=TRUE)

p2<-ggplot(m3a.pred.obs, aes(x=x, y=predicted))+
  geom_line(linewidth=1.5, color="darkgray")+geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=0.2)+
  stat_summary(fun=mean, fill="coral", color= "black", geom="point", size=5, shape=21)+
  geom_point(aes(x=x,y=prop.terrestrial), size=5,shape=21,stroke=0.8)+
  ylim(0,1)+
  theme(legend.title = element_blank(),panel.background = element_rect(fill = "white"),
        axis.title = element_text(size=16),
        axis.text = element_text(size=14),
        panel.border = element_rect(linetype = "solid", fill = NA))+
  ylab("Proportion of\nterrestrial species")+
  xlab("% Forest cover (500m)")


#########################
#    Wetland species    #
#########################

# Forest cover

m3b.pred<-ggpredict(m3b, terms = c("forest.500 [18.76,44.00,6.38,20.33,46.93,3.58,57.25,5.55,52.87,39.92,9.99,49.31,26.87,34.10,37.03,57.94,37.60,58.82,1.19,52.40,58.69,28.55]"))
write.csv(m3b.pred,"wet.forest.all.pred.csv",row.names=FALSE)
m3b.pred.obs<-read.csv("wet.forest.all.pred.csv", header=TRUE)

p3<-ggplot(m3b.pred.obs, aes(x=x, y=predicted))+
  geom_line(linewidth=1.5, color="darkgray")+geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=0.2)+
  stat_summary(fun=mean, fill="coral", color= "black", geom="point", size=5, shape=21)+
  geom_point(aes(x=x,y=prop.wetland), size=5,shape=21,stroke=0.8)+
  ylim(0,0.3)+
  theme(legend.title = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.title = element_text(size=16),
        axis.text = element_text(size=14),
        panel.border = element_rect(linetype = "solid", fill = NA))+
  ylab("Proportion of wetland \n species")+ xlab("% Forest Cover (500 m)")

#########################
#    Aquatic species    #
#########################

# Developed land

m4c.pred<-ggpredict(m4c, terms = c("developed.500 [46.26,6.33,72.81,7.39,9.08,17.68,24.04,81.25,33.09,22.95,49.12,7.20,20.18,10.68,19.82,9.15,28.32,13.98,92.20,21.65,12.42,1.50]"))
write.csv(m4c.pred,"aqu.dev.all.pred.csv",row.names=FALSE)
m4c.pred.obs<-read.csv("aqu.dev.all.pred.csv", header=TRUE)


p4<-ggplot(m4c.pred.obs, aes(x=x, y=predicted))+
  geom_line(linewidth=1.5,color="darkgray")+
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=0.2)+
  stat_summary(fun=mean, fill="coral", color= "black", geom="point", size=5, shape=21)+
  geom_point(aes(x=x,y=prop.aquatic), size=5,shape=21,stroke=0.8)+
  ylim(0,1)+
  theme(legend.title = element_blank(),panel.background = element_rect(fill = "white"),
        axis.title = element_text(size=16),
        axis.text = element_text(size=14),
        panel.border = element_rect(linetype = "solid", fill = NA))+
  ylab("Proportion of aquatic \n species")+
  xlab("% Developed land (500m)")


# Create the multipanel figure (Figure S3)
combinedPlot1 <- plot_grid(p1,p2,p3,p4, labels='auto', label_y= 1, label_x=0, align='hv', label_size=20, 
                           ncol=2, nrow = 2, rel_widths=1, rel_heights= 1)

png("Proportions regression models.png",width = 4000, height = 2200,
    units = "px", res=300)
combinedPlot1
dev.off()

