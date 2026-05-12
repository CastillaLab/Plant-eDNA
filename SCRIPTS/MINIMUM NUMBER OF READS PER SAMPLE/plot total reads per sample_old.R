library(ggplot2)
library(cowplot)
library(dplyr)

##Entire community

df1<-read.csv("FINAL_plants.median.clean.CM_0diff_NO_gen.rev_blast_28Dec21.csv", header=TRUE)
spinfo <- read.csv("spinfo.csv")

keep_spp<-spinfo$species
df1<-df1[,c(1:6,which(colnames(df1) %in% keep_spp))]



df1$total.reads = rowSums(df1[,7:227])

df1$lake_name <- gsub("_", " ", df1$lake_name)

df1$lake_name<-as.factor(df1$lake_name)
df1$sample<-as.factor(df1$sample)
df1$sample_type<-as.factor(df1$sample_type)

p1<-ggplot(df1,aes(x=total.reads))+
  geom_histogram(binwidth = 500, color="black",
                 fill="white")+
  theme_classic() + xlab("Total number of reads")+
  ylim(0,250)+
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size = 16))+
  geom_vline(xintercept = 1000, linetype="dashed",
             color="red", size=1)
  

p2<-ggplot(df1,aes(x=lake_name,y=total.reads))+
  geom_boxplot()+theme_classic()+ theme(axis.text.x = element_text(angle = 60,hjust = 1, size=8),
                                        axis.text.y = element_text(size=14),
                                        axis.title = element_text(size=16))+
  ylab("Total number \n of sequence reads")+
  xlab("Lake")+
  ylim(0,30000)

combinedPlot<-plot_grid(p1,p2,ncol=2,nrow=1,labels='auto',label_size = 20,
                        align='hv',rel_widths=1,rel_heights=1)


png("number of reads per sample.png",width = 4000, height = 1600, units = "px", res=300)
combinedPlot
dev.off()


