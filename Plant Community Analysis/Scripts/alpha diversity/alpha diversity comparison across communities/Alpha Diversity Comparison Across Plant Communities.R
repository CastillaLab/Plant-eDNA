################################################################################
# Script: Comparison of alpha diversity among plant community types
#
# Description:
# This script evaluates whether alpha diversity differs among
# terrestrial, wetland, and aquatic plant communities.
#
# Two alpha-diversity metrics are analyzed:
#   1. Species richness
#   2. Shannon diversity index
#
# The workflow:
#   1. Imports alpha-diversity data
#   2. Formats categorical variables
#   3. Explores distributions of diversity metrics
#   4. Fits generalized linear models (GLMs)
#   5. Compares fitted models against null models using
#      likelihood ratio tests
#   6. Generates plots
#
################################################################################


#######################################
# Load required packages
#######################################

# lmtest:
#   Provides statistical tests for model comparison.
#   In this workflow, lrtest() is used to compare
#   fitted GLMs against null models using
#   likelihood ratio tests.
library(lmtest)

# ggplot2:
#   Creates plots.
library(ggplot2)


################################################################################
# 1. IMPORT AND FORMAT DATA
################################################################################

#######################################
# Read alpha-diversity dataset
#######################################
df<-read.csv("alpha diversity lakes.csv", header=TRUE)


#######################################
# Convert categorical variables to factors
#######################################

# Plant community type
df$community<-as.factor(df$community)

# Lake identity
df$lake<-as.factor(df$lake)

#######################################
# Explore response variable distributions
#######################################

# Distribution of species richness
hist(df$richness)

# Distribution of Shannon diversity
hist(df$shannon)

################################################################################
# 2. SPECIES RICHNESS ANALYSIS
################################################################################

#######################################
# Fit generalized linear model
#######################################

# Model species richness as a function
# of plant community type.
#
# family = poisson is appropriate for
# count data such as species richness.
m1<-glm(richness~community,family=poisson,data=df)

# Display model summary
summary(m1)

#######################################
# Fit null model
#######################################

# Null model without community effects
m1.null<-glm(richness~1,family=poisson,data=df)

#######################################
# Compare models using likelihood ratio test
#######################################

# Test whether plant community type
# significantly improves model fit
lrtest(m1,m1.null)



################################################################################
# 3. VISUALIZATION OF SPECIES RICHNESS
################################################################################

#######################################
# Generate boxplot figure
#######################################
png("Richness among plant communities.png",width = 1500, height = 1000, units = "px", res=300)
ggplot(df,aes(x=community,y=richness,fill=community))+
  geom_boxplot(width=0.3, outlier.shape = NA)+
  geom_jitter(color="black", size=1, width=0.2, alpha=0.6)+
  theme_classic()+
  ylab("Species richness")+
  xlab("Plant community")+
  scale_fill_brewer(palette="Set2")+
  theme(axis.title.x = element_text(size=14,margin = margin(t = 15)),
        axis.title.y = element_text(size=14,margin = margin(r = 15)),
        axis.text = element_text(size=14),
        legend.position = "none")
dev.off()

################################################################################
# 4. SHANNON DIVERSITY ANALYSIS
################################################################################

#######################################
# Fit generalized linear model
#######################################

# Model Shannon diversity as a function
# of plant community type.
#
# Default Gaussian error structure is used
# because Shannon diversity is continuous.
m2<-glm(shannon~community,data=df)

# Display model summary
summary(m2)

#######################################
# Fit null model
#######################################

# Null model without community effects
m2.null<-glm(shannon~1,data=df)

#######################################
# Compare models using likelihood ratio test
#######################################

# Test whether plant community type
# significantly improves model fit
lrtest(m2,m2.null)

