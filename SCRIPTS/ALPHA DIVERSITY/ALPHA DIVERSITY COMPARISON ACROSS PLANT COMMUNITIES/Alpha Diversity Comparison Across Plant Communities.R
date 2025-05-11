library(lmtest)
library(ggplot2)

df<-read.csv("alpha diversity lakes.csv", header=TRUE)

df$community<-as.factor(df$community)
df$lake<-as.factor(df$lake)

hist(df$richness)
hist(df$shannon)

#############################
# Species richness analysis #
#############################

# Create a model comparing species richness among terrestrial, 
# wetland, and aquatic communities
m1<-glm(richness~community,family=poisson,data=df)
summary(m1)

# Create the null model
m1.null<-glm(richness~1,family=poisson,data=df)

# Compare the model with community type with the null model using the 
# Likelihood Ratio test
lrtest(m1,m1.null)

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

#############################
#     Shannon Diversity     #
#############################

# Create a model comparing Shannon diversity among terrestrial, 
# wetland, and aquatic communities
m2<-glm(shannon~community,data=df)
summary(m2)

# Create the null model
m2.null<-glm(shannon~1,data=df)

# Compare the model with community type with the null model using the 
# Likelihood Ratio test
lrtest(m2,m2.null)

