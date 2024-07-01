library(ggplot2)
library(cowplot)

# Plot with proportions of reads of terrestrial, wetland, and aquatic plant
# species

df<-read.csv("reads.richness.groups.csv",header=TRUE)

p1<-ggplot(df,aes(x=terrestrial.spp, y=aquatic.spp))+
  geom_smooth(method = "lm", se=FALSE, col = "black", size=0.8)+
  geom_point(shape=21,size=3)+theme_bw()+
  xlab("Richness Terrestrial Species")+
  ylab("Richness Aquatic \nSpecies")+
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=14))


p2<-ggplot(df,aes(x=wetland.spp, y=aquatic.spp))+
  geom_smooth(method = "lm", se=FALSE, col = "black", size=0.8)+
  geom_point(shape=21,size=3)+theme_bw()+
  xlab("Richness Wetland Species")+
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_blank(),
        axis.text = element_text(size=14))

p3<-ggplot(df,aes(x=terrestrial.reads, y=aquatic.reads))+
  geom_smooth(method = "lm", se=FALSE, col = "black", size=0.8)+
  geom_point(shape=21,size=3)+theme_bw()+
  xlab("Read counts for \nterrestrial spp")+
  ylab("Read counts for \naquatic spp")+
  scale_x_continuous(labels = scales::comma)+
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=14))

p4<-ggplot(df,aes(x=wetland.reads, y=aquatic.reads))+
  geom_smooth(method = "lm", se=FALSE, col = "black", size=0.8)+
  geom_point(shape=21,size=3)+theme_bw()+
  xlab("Read counts for \nwetland spp")+
  ylab("Read counts for aquatic spp")+
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_blank(),
        axis.text = element_text(size=14))

p5<-ggplot(df,aes(x=prop.aquatic, y=prop.terrestrial))+
  geom_smooth(method = "lm", se=FALSE, col = "black", size=0.8)+
  geom_point(shape=21,size=3)+theme_bw()+
  xlab("Proportion of reads for \n terrestrial spp")+
  ylab("Proportion of reads for \n aquatic spp")+
  xlim(0,1)+ylim(0,1)+
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=14))


p6<-ggplot(df,aes(x=prop.aquatic, y=prop.wetland))+
  geom_smooth(method = "lm", se=FALSE, col = "black", size=0.8)+
  geom_point(shape=21,size=3)+theme_bw()+
  xlab("Proportion of reads for \n wetland spp")+
  ylab("Proportion of reads for \n aquatic spp")+
  xlim(0,1)+ylim(0,0.2)+
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_blank(),
        axis.text = element_text(size=14))

combinedPlot <- plot_grid(p1, p2, p3, p4,p5, p6,labels='auto',label_y= 1.02,
                          label_x=0.06,label_size=20,ncol=2, nrow=3, align='hv',
                          rel_widths=1,rel_heights= 1)
                         

png("Correlations reads and richness.png",width = 3500, height = 3000, units = "px", res=300)
combinedPlot
dev.off()


# Correlations

cor.test(df$aquatic.spp,df$terrestrial.spp)
cor.test(df$aquatic.spp,df$wetland.spp)

cor.test(df$aquatic.reads,df$terrestrial.reads)
cor.test(df$aquatic.reads,df$wetland.reads)

cor.test(df$prop.aquatic,df$prop.terrestrial, method = "spearman")
cor.test(df$prop.aquatic,df$prop.wetland, method = "spearman")

