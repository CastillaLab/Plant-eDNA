library(ggplot2)
library(cowplot)

# Plot with proportions of reads of terrestrial, wetland, and aquatic plant
# species

df<-read.csv("prop.reads.richness.groups.csv",header=TRUE)

df$lake.name<-as.factor(df$lake.name)
df$habitat<-as.factor(df$habitat)

df$lake.name <- factor(df$lake.name, levels = c("Long","Fourth","Austin","Dumont","Cass",
                                                  "Thompson","Holloway.Reservoir","Haithco","Pickerel",
                                                  "Kimball","Pentwater","George","Houghton",
                                                  "Five.Channels","Higgins","Torch","Walloon",
                                                  "Ocqueoc","Mullett","Wycamp","Brevort","Big.Manistique"))

p1<-ggplot(df, aes(fill=habitat, y=proportion, x=lake.name)) + 
  geom_bar(position="stack", stat="identity") +
  theme_classic()+theme(axis.text.x = element_text(angle = 60,hjust = 1, size=8),
                        axis.text.y = element_text(size=12),
                        legend.text=element_text(size=rel(1)),
                        legend.title = element_text(size=rel(1.25)),
                        axis.title.x = element_blank(),
                        axis.title.y = element_text(size=18,
                                                    margin = margin(t=0,r=10,b=0,l=0)))+
  scale_fill_brewer(palette="Set2")+ylab("Proportion of reads")


# Plot species richness


p2<-ggplot(df, aes(fill=habitat, y=species.richness, x=lake.name)) + 
  geom_bar(position="stack", stat="identity") +
  theme_classic()+theme(axis.text.x = element_text(angle = 60,hjust = 1, size=8),
                        axis.text.y = element_text(size=12),
                        axis.title.x = element_blank(),
                        legend.position = "none",
                        axis.title.y = element_text(size=18,
                                                    margin = margin(t=0,r=10,b=0,l=0)))+
  scale_fill_brewer(palette="Set2")+ylab("Species richness")+ylim(0,100)


combinedPlot <- plot_grid(p2, p1, ncol=2, labels='auto',label_y= 1,
                          label_x=0,label_size=20,align='hv', rel_widths=1,
                          rel_heights= 1)

png("Proportions per lake.png",width = 3500, height = 1400, units = "px", res=300)
combinedPlot
dev.off()



