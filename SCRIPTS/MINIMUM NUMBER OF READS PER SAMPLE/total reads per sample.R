library(ggplot2)
library(cowplot)

##Entire community

df1<-read.csv("eDNA plant community_entire.csv", header=TRUE)

df1$total.reads = rowSums(df1[,4:215])


df1$lake.name<-as.factor(df1$lake.name)
df1$sample<-as.factor(df1$sample)
df1$sample.type<-as.factor(df1$sample.type)

p1<-ggplot(df1,aes(x=total.reads))+
  geom_histogram(binwidth = 500, color="black",
                 fill="white")+
  theme_classic() + xlab("Total number of reads")+
  ylim(0,250)+
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size = 16))+
  geom_vline(xintercept = 1000, linetype="dashed",
             color="red", size=1)
  

p2<-ggplot(df1,aes(x=lake.name,y=total.reads))+
  geom_boxplot()+theme_classic()+ theme(axis.text.x = element_text(angle = 60,hjust = 1, size=8),
                                        axis.text.y = element_text(size=14),
                                        axis.title = element_text(size=16))+
  ylab("Total number \n of sequence reads")+
  xlab("Lake")+
  ylim(0,30000)

combinedPlot<-plot_grid(p1,p2,ncol=2,nrow=1,labels='auto',label_size = 20,
                        align='hv',rel_widths=1,rel_heights=1)


png("number of reads per sample.png",width = 3000, height = 1600, units = "px", res=300)
combinedPlot
dev.off()


