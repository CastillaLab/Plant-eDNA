#Creating the species richness counts to use in heatmaps
setwd("~/Downloads")
cm <- read.csv("FINAL_plants.median.clean.CM_0diff_NO_gen.rev_blast_28Dec21_metadata.csv", header = TRUE)
spinfo <- read.csv("spinfo.csv", header = TRUE)

out <- cm[,c(1,2,3,4,16,17)]

spcntmat <- cm[,-c(1:17)]
spcntmat <- spcntmat[,-grep("unclassified",colnames(spcntmat))]

tot <- rowSums(spcntmat > 0)
aq <- rowSums(spcntmat[,which(colnames(spcntmat) %in% spinfo$species[spinfo$habitat == "aquatic"])] > 0)
wet <- rowSums(spcntmat[,which(colnames(spcntmat) %in% spinfo$species[spinfo$habitat == "wetland"])] > 0)
terr <- rowSums(spcntmat[,which(colnames(spcntmat) %in% spinfo$species[spinfo$habitat == "terrestrial"])] > 0)

out$Total <- tot
out$Aquatic <- aq
out$Wetland <- wet
out$Terrestrial <- terr

out

write.csv(out, "SRcounts_byHab.csv", quote = FALSE, row.names = FALSE)