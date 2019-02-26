# install.packages("ggbiplot")
# install.packages("heatmap3")
# install.packages("phyloseq")
# install.packages("vegan")
# install.packages("DESeq2")
# install.packages("phyloseq")
# 
# source("https://bioconductor.org/biocLite.R")
# biocLite("phyloseq")
# biocLite("ade4")

require("heatmap3")
library("data.table")
require(plyr)
require(vegan)
library(reshape2)
library("phyloseq")
library("DESeq2")
library("genefilter")
library("ape")
library("phytools")
library("vegan")
library("ggplot2")
library("plot3D")
library("calibrate")
library("ggbiplot")
library(reshape)
# library(ggplot2)
# library(limma)
# require(cowplot)
# library(gridGraphics)
dev.off()
rm(list=ls())


setwd("/Users/ahmedmetwally/Dropbox/LungTransplantWASHU/CF_Paper_v4.0")
setwd("D:/Dropbox/LungTransplantWASHU/CF_Paper_v4.0")



#####################################
#### Summarize MetaData #############
#####################################
meta = read.csv(file="metadata/CF_12_v2.4.csv", header = TRUE)


# ###### Add columns to the metaData
# ## Create sample ID
# meta$SUBJECT.WUTX = as.factor(meta$SUBJECT.WUTX)
# meta$SUBJECT.WUTX = gsub("SUBJ", "S", meta$SUBJECT.WUTX)
# meta$MiSeqRunNum = as.factor(meta$MiSeqRunNum)
# meta$id = paste(meta$SUBJECT.WUTX, meta$TimePoint, sep = ".")
# 
# ### Binning
# meta$bin = 0
# meta$bin[meta$TimeAfterTransplantation<50] = 50
# meta$bin[abs(meta$TimeAfterTransplantation - 100) < 50] = 100
# meta$bin[abs(meta$TimeAfterTransplantation - 200) < 50] = 200
# meta$bin[abs(meta$TimeAfterTransplantation - 300) < 50] = 300
# meta$bin[abs(meta$TimeAfterTransplantation - 400) < 50] = 400
# meta$bin[abs(meta$TimeAfterTransplantation - 500) < 50 | meta$TimeAfterTransplantation>500] = 500
# 
# ## Save the updated MetaFile
# #write.csv(meta, file = "CF_12_v1.1.csv", row.names = FALSE)



cbbPalette <- c("#999999", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000", "#E69F00", "#56B4E9", "#009E73", "darkorchid",
                "mistyrose1", "sienna1")
cbPalette <- c("#D55E00", "#009E73", "#56B4E9")

#####################################
#### Summarize Sequencing data ######
#####################################
meta$id_order = paste(meta$AnnotatedID, meta$TimePoint, sep = ".")
###################
### Raw Reads
###################
jpeg("LgTxCF12_RawSeq_small.jpg", res = 300, height = 15, width = 30, units = 'cm')
p = ggplot(meta, aes(x = id_order, y = RawReads, fill = AnnotatedID)) + theme_bw() 
p + geom_bar(stat="identity", position="stack") + theme(legend.position="bottom") + ggtitle("# raw reads per sample") + 
  labs(y = "# reads", x = "Sample") +
  scale_fill_manual(values = cbbPalette) +
  theme(axis.text.x = element_text(colour="black", size=10, angle=90, hjust=0.5, vjust=0, face="bold"),
        axis.text.y = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="bold"),
        axis.title.x = element_text(colour="black", size=15, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=15, angle=90, hjust=.5, vjust=.5, face="bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.text=element_text(size=15, face="plain"), legend.title = element_blank())
dev.off()



###################
### Screened Reads
###################
jpeg("LgTxCF12_Screened.jpg", res = 300, height = 15, width = 30, units = 'cm')
p = ggplot(meta, aes(x = id_order, y = ScreenedReads, fill = AnnotatedID)) + theme_bw() 
p + geom_bar(stat="identity", position="stack") + theme(legend.position="bottom") + ggtitle("# raw screened reads per sample") + 
  labs(y = "# reads", x = "Sample") +
  scale_fill_manual(values = cbbPalette) +
  theme(axis.text.x = element_text(colour="black", size=10, angle=90, hjust=0.5, vjust=0, face="bold"),
        axis.text.y = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="bold"),
        axis.title.x = element_text(colour="black", size=15, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=15, angle=90, hjust=.5, vjust=.5, face="bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.text=element_text(size=15, face="plain"), legend.title = element_blank())
dev.off()



#####################################
### Combined raw and screened Reads
####################################
### TODO: Do it with Violin plot for each subject

df = meta[, c("id", "RawReads", "ScreenedReads", "NumContigs")]
df.melt = melt(df, id.vars = "id")

jpeg("LgTxCF12_Filtered_Screened.jpg", res = 300, height = 20, width = 30, units = 'cm')
p = ggplot(df.melt, aes(x = id, y = value)) + theme_bw() 
p + geom_bar(aes(fill = variable), stat="identity", position="dodge") + theme(legend.position="bottom") + ggtitle("Screened read count") + 
  labs(y = "# reads", x = "Sample") +
  scale_fill_manual(values = cbPalette) +
  theme(axis.text.x = element_text(colour="black", size=10, angle=90, hjust=0.5, vjust=0, face="bold"),
        axis.text.y = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="bold"),
        axis.title.x = element_text(colour="black", size=15, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=15, angle=90, hjust=.5, vjust=.5, face="bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.text=element_text(size=15, face="plain"), legend.title = element_blank())
dev.off()



#####################################
####### Contigs number
####################################
jpeg("LgTxCF12_ContigsNumber.jpg", res = 300, height = 20, width = 30, units = 'cm')
p = ggplot(meta, aes(x = id, y = NumContigs, fill = SUBJECT.WUTX)) + theme_bw() 
p + geom_bar(stat="identity", position="stack") + theme(legend.position="bottom") + ggtitle("# Contigs") + 
  labs(y = "# contigs", x = "Sample") +
  scale_fill_manual(values = cbbPalette) +
  theme(axis.text.x = element_text(colour="black", size=10, angle=90, hjust=0.5, vjust=0, face="bold"),
        axis.text.y = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="bold"),
        axis.title.x = element_text(colour="black", size=15, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=15, angle=90, hjust=.5, vjust=.5, face="bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.text=element_text(size=15, face="plain"), legend.title = element_blank())
dev.off()


###########################
###### Contigs N50 ########
###########################
jpeg("LgTxCF12_N50.jpg", res = 300, height = 20, width = 30, units = 'cm')
p = ggplot(meta, aes(x = id, y = N50, fill = SUBJECT.WUTX)) + theme_bw() 
p + geom_bar(stat="identity", position="stack") + theme(legend.position="bottom") + ggtitle("Contigs' N50") + 
  labs(y = "N50 Length", x = "Sample") +
  scale_fill_manual(values = cbbPalette) +
  theme(axis.text.x = element_text(colour="black", size=10, angle=90, hjust=0.5, vjust=0, face="bold"),
        axis.text.y = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="bold"),
        axis.title.x = element_text(colour="black", size=15, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=15, angle=90, hjust=.5, vjust=.5, face="bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.text=element_text(size=15, face="plain"), legend.title = element_blank())
dev.off()


################################
###### Longest Contigs Length ##
################################
jpeg("LgTxCF12_longest_contig.jpg", res = 300, height = 20, width = 30, units = 'cm')
p = ggplot(meta, aes(x = id, y = LongestContig, fill = SUBJECT.WUTX)) + theme_bw() 
p + geom_bar(stat="identity", position="stack") + theme(legend.position="bottom") + ggtitle("Longest contig length") + 
  labs(y = "Contig length (bp)", x = "Sample") +
  scale_fill_manual(values = cbbPalette) +
  theme(axis.text.x = element_text(colour="black", size=10, angle=90, hjust=0.5, vjust=0, face="bold"),
        axis.text.y = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="bold"),
        axis.title.x = element_text(colour="black", size=15, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=15, angle=90, hjust=.5, vjust=.5, face="bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.text=element_text(size=15, face="plain"), legend.title = element_blank())
dev.off()


###################################
#### Visualize timepoints dist
####################################
TimePoints = data.frame(PatientID = as.factor(paste("",meta$AnnotatedID, sep="")), DaysAfterTransplantation = meta$TimeAfterTransplantation, Group= meta$OutcomeMod)
TimePoints = TimePoints[order(TimePoints$Group, TimePoints$PatientID),]
TimePoints = data.frame(PatientID = factor(TimePoints$PatientID, levels = TimePoints$PatientID), DaysAfterTransplantation = TimePoints$DaysAfterTransplantation, Group= as.factor(TimePoints$Group))


jpeg("LgTxCF12_timepoints_distribution.jpg", res = 300, height = 10, width = 20, units = 'cm')
ggplot(TimePoints, aes(x = PatientID, y = DaysAfterTransplantation)) +
  theme_bw() + 
  scale_colour_manual(values=c( "red", "midnightblue"),  breaks=c("BO", "non-BO")) +
  # theme(axis.text.x=element_text(angle=45, hjust=1)) +
  scale_x_discrete(name ="Subject")+ #coord_fixed(ratio = 0.002)+
  geom_point(aes(colour = Group), stat='identity', size = 2, show.legend = TRUE) +
  theme(legend.position="top")+ scale_shape_manual(values=c(1,6)) +
  scale_y_continuous(name= "Days after Transplantation", breaks = round(seq(0, max(TimePoints$DaysAfterTransplantation) + 20, by = 100),1)) +
  theme(axis.text.x = element_text(colour="black",size=12,angle=45,hjust=1,vjust=1,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=1,vjust=0.5,face="plain"),
        axis.title.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=0.5,face="bold"),
        axis.title.y = element_text(colour="black",size=12,angle=90,hjust=.5,vjust=.5,face="bold"),
        legend.text=element_text(size=15, face="bold"), legend.title = element_blank()) +
  coord_flip()
# h = recordPlot(load="gridGraphics", attach=NULL)
dev.off()



#######################################
#### Visualize Binned timepoints dist
#######################################
TimePoints = data.frame(PatientID = as.factor(meta$AnnotatedID), DaysAfterTransplantation = meta$bin, Group= as.factor(meta$Outcome))
TimePoints = TimePoints[order(TimePoints$Group, TimePoints$PatientID),]
TimePoints = data.frame(PatientID = factor(TimePoints$PatientID, levels = TimePoints$PatientID), DaysAfterTransplantation = TimePoints$DaysAfterTransplantation, Group= as.factor(TimePoints$Group))

jpeg("LgTxCF12_timepoints_distribution_binned.jpg", res = 300, height = 10, width = 10, units = 'cm')
ggplot(TimePoints, aes(x = PatientID , y = DaysAfterTransplantation)) +
  theme_bw() + 
  scale_colour_manual(values=c( "red", "midnightblue"),  breaks=c("BO", "non-BO")) +
  # theme(axis.text.x=element_text(angle=45, hjust=1)) +
  scale_x_discrete(name ="Subject")+ #coord_fixed(ratio = 0.002)+
  geom_point(aes(colour = Group), stat='identity', size = 2, show.legend = TRUE, alpha = 0.6) +
  theme(legend.position="top")+ scale_shape_manual(values=c(1,6)) +
  scale_y_continuous(name= "Days after Transplantation", breaks = round(seq(0, max(TimePoints$DaysAfterTransplantation) + 20, by = 100),1)) +
  theme(axis.text.x = element_text(colour="black",size=10,angle=45,hjust=1,vjust=1,face="plain"),
        axis.text.y = element_text(colour="black",size=10,angle=0,hjust=1,vjust=0.5,face="plain"),
        axis.title.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=0.5,face="bold"),
        axis.title.y = element_text(colour="black",size=12,angle=90,hjust=.5,vjust=.5,face="bold"),
        legend.text=element_text(size=15, face="bold"), legend.title = element_blank()) +
  coord_flip()
# h = recordPlot(load="gridGraphics", attach=NULL)
dev.off()

  


#######################################################################
#### Boxplot for FEV and FVC per subject and per outcome
#######################################################################
#par(mfrow=c(1,3))

PFT = data.frame(PatientID = as.factor(meta$SUBJECT.WUTX), FVC = meta$FVC_Per_Predicted, FEV = meta$FEV1_Per_Predicted, Group= as.factor(meta$Outcome))
PFT = PFT[order(PFT$Group, PFT$PatientID),]
PFT = data.frame(PatientID = factor(PFT$PatientID, levels = PFT$PatientID), FVC = PFT$FVC, FEV = PFT$FEV, Group= as.factor(PFT$Group))


######################
###### FEV ###########
######################
jpeg("LgTxCF12_FEV1.jpg", res = 300, height = 10, width = 15, units = 'cm')
ggplot(PFT, aes(x = PatientID, y = FEV, fill = factor(Group))) +
  #geom_violin(aes(fill = Outcome), trim = FALSE) + 
  geom_boxplot() +
  geom_dotplot(binaxis = 'y', stackdir='center', dotsize = 0.5, position = position_dodge(0.2), color = "black", fill = "white") + 
  #geom_dotplot(binaxis = 'y', stackdir='center', dotsize = 0.5, position = position_dodge(0.2), aes(fill = TimeAfterTransplantation)) + 
  ggtitle("FEV1 (Percent of Perdicted)") + 
  scale_fill_manual(values=c( "red", "midnightblue"),  breaks=c("bos", "non-bos")) + 
  #scale_colour_manual(values=c( "red", "midnightblue"),  breaks=c("bos", "non-bos")) +
  #scale_x_discrete(name ="Subject") +
  #scale_y_discrete(name ="FVC Predicted (%)") +
  #scale_fill_manual(values = cbPalette) +
  #scale_color_manual(values = cbPalette) +
  #scale_fill_brewer(palette = "Set3")
  labs(y = "FEV1 % of Predicted", x = "Subject") + 
  theme_bw() + #stat_summary(fun.y=mean, geom="point", size=2, color="black") + geom_boxplot(width=0.1) +
  #scale_y_continuous(breaks = c(0, 1)) + 
  theme(axis.text.x = element_text(colour="black", size=10, angle=45, hjust=0, vjust=0, face="plain"),
        axis.text.y = element_text(colour="black", size=10, angle=0, hjust=0.5, vjust=0.5, face="plain"),
        axis.title.x = element_text(colour="black", size=12, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=12, angle=90, hjust=.5, vjust=.5, face="bold"),
        #legend.text = element_text(size=15, face="plain"), 
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position="top",
        strip.text = element_text(size=12, face = "bold"),
        legend.text=element_text(size=14, face="bold"), 
        legend.title = element_blank())
dev.off()


######################
###### FVC  ##########
######################
jpeg("LgTxCF12_FVC.jpg", res = 300, height = 10, width = 15, units = 'cm')
ggplot(PFT, aes(x = PatientID, y = FVC, fill = factor(Group))) +
  #geom_violin(aes(fill = Outcome), trim = FALSE) + 
  geom_boxplot() +
  geom_dotplot(binaxis = 'y', stackdir='center', dotsize = 0.5, position = position_dodge(0.2), color = "black", fill = "white") + 
  #geom_dotplot(binaxis = 'y', stackdir='center', dotsize = 0.5, position = position_dodge(0.2), aes(fill = TimeAfterTransplantation)) + 
  ggtitle("FVC (Percent of Perdicted)") + 
  scale_fill_manual(values=c( "red", "midnightblue"),  breaks=c("bos", "non-bos")) + 
  #scale_colour_manual(values=c( "red", "midnightblue"),  breaks=c("bos", "non-bos")) +
  #scale_x_discrete(name ="Subject") +
  #scale_y_discrete(name ="FVC Predicted (%)") +
  #scale_fill_manual(values = cbPalette) +
  #scale_color_manual(values = cbPalette) +
  #scale_fill_brewer(palette = "Set3")
  labs(y = "FVC % of Predicted", x = "Subject") + 
  theme_bw() + #stat_summary(fun.y=mean, geom="point", size=2, color="black") + geom_boxplot(width=0.1) +
  #scale_y_continuous(breaks = c(0, 1)) + 
  theme(axis.text.x = element_text(colour="black", size=10, angle=45, hjust=0, vjust=0, face="plain"),
        axis.text.y = element_text(colour="black", size=10, angle=0, hjust=0.5, vjust=0.5, face="plain"),
        axis.title.x = element_text(colour="black", size=12, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=12, angle=90, hjust=.5, vjust=.5, face="bold"),
        #legend.text=element_text(size=15, face="plain"), 
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position="top",
        strip.text = element_text(size=12, face = "bold"),
        legend.text=element_text(size=14, face="bold"), 
        legend.title = element_blank())
dev.off()



#########################
#### Overall FEV1
########################
jpeg("LgTxCF12_FEV1_overall.jpg", res = 300, height = 12, width = 12, units = 'cm')
ggplot(meta, aes(x = factor(Outcome), y = FEV1_PercPredicted, fill = factor(Outcome))) +
  geom_violin(aes(fill = Outcome), trim = FALSE) + 
  #geom_boxplot() +
  geom_dotplot(binaxis = 'y', stackdir='center', dotsize = 0.7, position = position_dodge(.2), color = "black", fill = "white") + 
  #geom_dotplot(binaxis = 'y', stackdir='center', dotsize = 0.5, position = position_dodge(0.2), aes(fill = TimeAfterTransplantation)) + 
  ggtitle("FEV1 (Percent of Perdicted) p-value=0.23") + 
  scale_fill_manual(values=c( "red", "midnightblue"),  breaks=c("bos", "non-bos")) + 
  #scale_colour_manual(values=c( "red", "midnightblue"),  breaks=c("bos", "non-bos")) +
  #scale_x_discrete(name ="Subject") +
  #scale_y_discrete(name ="FVC Predicted (%)") +
  #scale_fill_manual(values = cbPalette) +
  #scale_color_manual(values = cbPalette) +
  #scale_fill_brewer(palette = "Set3")
  labs(y = "FEV1 % of Predicted", x = "Subject") + 
  theme_bw() + #stat_summary(fun.y=mean, geom="point", size=2, color="black") + geom_boxplot(width=0.1) +
  #scale_y_continuous(breaks = c(0, 1)) + 
  theme(axis.text.x = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0, face="bold"),
        axis.text.y = element_text(colour="black", size=10, angle=0, hjust=0.5, vjust=0.5, face="plain"),
        axis.title.x = element_blank(), #element_text(colour="black", size=12, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=12, angle=90, hjust=.5, vjust=.5, face="bold"),
        #legend.text=element_text(size=15, face="plain"), 
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position="top",
        strip.text = element_text(size=12, face = "bold"),
        legend.text=element_text(size=14, face="bold"), 
        legend.title = element_blank())
dev.off()



#########################
#### Overall FVC
########################
jpeg("LgTxCF12_FVC_overall.jpg", res = 300, height = 12, width = 12, units = 'cm')
ggplot(meta, aes(x = factor(Outcome), y = meta$FVC_PercPredicted, fill = factor(Outcome))) +
  geom_violin(aes(fill = Outcome), trim = FALSE) + 
  #geom_boxplot() +
  geom_dotplot(binaxis = 'y', stackdir='center', dotsize = 0.7, position = position_dodge(.2), color = "black", fill = "white") + 
  #geom_dotplot(binaxis = 'y', stackdir='center', dotsize = 0.5, position = position_dodge(0.2), aes(fill = TimeAfterTransplantation)) + 
  ggtitle("FVC (Percent of Perdicted) p-value=0.16") + 
  scale_fill_manual(values=c( "red", "midnightblue"),  breaks=c("bos", "non-bos")) + 
  #scale_colour_manual(values=c( "red", "midnightblue"),  breaks=c("bos", "non-bos")) +
  #scale_x_discrete(name ="Subject") +
  #scale_y_discrete(name ="FVC Predicted (%)") +
  #scale_fill_manual(values = cbPalette) +
  #scale_color_manual(values = cbPalette) +
  #scale_fill_brewer(palette = "Set3")
  labs(y = "FVC % of Predicted", x = "Subject") + 
  theme_bw() + #stat_summary(fun.y=mean, geom="point", size=2, color="black") + geom_boxplot(width=0.1) +
  #scale_y_continuous(breaks = c(0, 1)) + 
  theme(axis.text.x = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0, face="bold"),
        axis.text.y = element_text(colour="black", size=10, angle=0, hjust=0.5, vjust=0.5, face="plain"),
        axis.title.x = element_blank(), #element_text(colour="black", size=12, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=12, angle=90, hjust=.5, vjust=.5, face="bold"),
        #legend.text=element_text(size=15, face="plain"), 
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position="top",
        strip.text = element_text(size=12, face = "bold"),
        legend.text=element_text(size=14, face="bold"), 
        legend.title = element_blank())
dev.off()

#######################################
######## P-value for FEV and FVC:  ####
#######################################
wilcox.test(FVC_Per_Predicted ~ Outcome, data = meta)
#data:  FVC_Per_Predicted by Outcome
#W = 457.5, p-value = 0.09012

wilcox.test(FEV1_Per_Predicted ~ Outcome, data = meta)
#data:  FEV1_Per_Predicted by Outcome
#W = 479, p-value = 0.1493




#################################
### Extract metadata for Table
#################################
metaDemog = data.frame(meta$SUBJECT.WUTX, meta$NumPoints, meta$RecipientAge, meta$RecipientGender, 
                       meta$FEV1_beforeTx, meta$RightLungPerfusionFindings, meta$Outcome, meta$RecipientWeightKG, 
                       meta$Genotyping, meta$GenotypingClass, meta$BMI_kg_m2)
metaDemog_uniq = unique(metaDemog)



### GenoTyping
# metaDemog_uniq[which(metaDemog_uniq$meta.Genotyping == "DF508/DF508"), "meta.Genotyping"] = factor("homo")
# 
# 
# metaDemog_uniq[metaDemog_uniq$meta.Genotyping == "DF508/DF508", "meta.Genotyping"] = factor("homo")

##  Gender
table(metaDemog_uniq[metaDemog_uniq$meta.Outcome == "bos",]$meta.RecipientGender)
table(metaDemog_uniq[metaDemog_uniq$meta.Outcome == "non-bos",]$meta.RecipientGender)



## Age
summary(metaDemog_uniq$meta.RecipientAge)
sd(metaDemog_uniq$meta.RecipientAge)

summary(metaDemog_uniq[metaDemog_uniq$meta.Outcome == "bos",]$meta.RecipientAge)
sd(metaDemog_uniq$meta.RecipientAge)

summary(metaDemog_uniq[metaDemog_uniq$meta.Outcome == "non-bos",]$meta.RecipientAge)
sd(metaDemog_uniq$meta.RecipientAge)


wilcox.test(meta.RecipientAge ~ meta.Outcome, data = metaDemog_uniq)

## Weight
mean(metaDemog_uniq$meta.RecipientWeightKG)
sd(metaDemog_uniq$meta.RecipientWeightKG)
summary(metaDemog_uniq$meta.RecipientWeightKG)

mean(metaDemog_uniq[metaDemog_uniq$meta.Outcome == "bos",]$meta.RecipientWeightKG)
sd(metaDemog_uniq[metaDemog_uniq$meta.Outcome == "bos",]$meta.RecipientWeightKG)
summary(metaDemog_uniq[metaDemog_uniq$meta.Outcome == "bos",]$meta.RecipientWeightKG)

mean(metaDemog_uniq[metaDemog_uniq$meta.Outcome == "non-bos",]$meta.RecipientWeightKG)
sd(metaDemog_uniq[metaDemog_uniq$meta.Outcome == "non-bos",]$meta.RecipientWeightKG)
summary(metaDemog_uniq[metaDemog_uniq$meta.Outcome == "non-bos",]$meta.RecipientWeightKG)

wilcox.test(meta.RecipientWeightKG ~ meta.Outcome, data = metaDemog_uniq)



## BMI
mean(metaDemog_uniq$meta.BMI_kg_m2)
sd(metaDemog_uniq$meta.BMI_kg_m2)
summary(metaDemog_uniq$meta.BMI_kg_m2)

mean(metaDemog_uniq[metaDemog_uniq$meta.Outcome == "bos",]$meta.BMI_kg_m2)
sd(metaDemog_uniq[metaDemog_uniq$meta.Outcome == "bos",]$meta.BMI_kg_m2)
summary(metaDemog_uniq[metaDemog_uniq$meta.Outcome == "bos",]$meta.BMI_kg_m2)

mean(metaDemog_uniq[metaDemog_uniq$meta.Outcome == "non-bos",]$meta.BMI_kg_m2)
sd(metaDemog_uniq[metaDemog_uniq$meta.Outcome == "non-bos",]$meta.BMI_kg_m2)
summary(metaDemog_uniq[metaDemog_uniq$meta.Outcome == "non-bos",]$meta.BMI_kg_m2)

wilcox.test(meta.BMI_kg_m2 ~ meta.Outcome, data = metaDemog_uniq)




#### TimePoints 
summary(metaDemog_uniq$meta.NumPoints)
sd(metaDemog_uniq$meta.NumPoints)

summary(metaDemog_uniq[metaDemog_uniq$meta.Outcome == "bos",]$meta.NumPoints)
sd(metaDemog_uniq$meta.NumPoints)

summary(metaDemog_uniq[metaDemog_uniq$meta.Outcome == "non-bos",]$meta.NumPoints)
sd(metaDemog_uniq$meta.NumPoints)

wilcox.test(meta.NumPoints ~ meta.Outcome, data = metaDemog_uniq)


## FEV1_beforeTx
mean(metaDemog_uniq$meta.FEV1_beforeTx)
sd(metaDemog_uniq$meta.FEV1_beforeTx)
summary(metaDemog_uniq$meta.FEV1_beforeTx)

mean(metaDemog_uniq[metaDemog_uniq$meta.Outcome == "bos",]$meta.FEV1_beforeTx)
sd(metaDemog_uniq[metaDemog_uniq$meta.Outcome == "bos",]$meta.FEV1_beforeTx)
summary(metaDemog_uniq[metaDemog_uniq$meta.Outcome == "bos",]$meta.FEV1_beforeTx)

mean(metaDemog_uniq[metaDemog_uniq$meta.Outcome == "non-bos",]$meta.FEV1_beforeTx)
sd(metaDemog_uniq[metaDemog_uniq$meta.Outcome == "non-bos",]$meta.FEV1_beforeTx)
summary(metaDemog_uniq[metaDemog_uniq$meta.Outcome == "non-bos",]$meta.FEV1_beforeTx)

wilcox.test(meta.FEV1_beforeTx ~ meta.Outcome, data = metaDemog_uniq)


### Right Lung Perfusion
mean(metaDemog_uniq$meta.RightLungPerfusionFindings)
sd(metaDemog_uniq$meta.RightLungPerfusionFindings)
summary(metaDemog_uniq$meta.RightLungPerfusionFindings)

mean(metaDemog_uniq[metaDemog_uniq$meta.Outcome == "bos",]$meta.RightLungPerfusionFindings)
sd(metaDemog_uniq[metaDemog_uniq$meta.Outcome == "bos",]$meta.RightLungPerfusionFindings)
summary(metaDemog_uniq[metaDemog_uniq$meta.Outcome == "bos",]$meta.RightLungPerfusionFindings)

mean(metaDemog_uniq[metaDemog_uniq$meta.Outcome == "non-bos",]$meta.RightLungPerfusionFindings)
sd(metaDemog_uniq[metaDemog_uniq$meta.Outcome == "non-bos",]$meta.RightLungPerfusionFindings)
summary(metaDemog_uniq[metaDemog_uniq$meta.Outcome == "non-bos",]$meta.RightLungPerfusionFindings)

wilcox.test(meta.RightLungPerfusionFindings ~ meta.Outcome, data = metaDemog_uniq)






##### Categorical Variable
male = metaDemog_uniq[metaDemog_uniq$meta.RecipientGender == "Male",]
tbl = table(male$meta.RecipientGender, male$meta.Outcome)
chisq.test(tbl)

tbl = table(male$meta.RecipientGender, male$meta.Outcome)
chisq.test(tbl)

tbl = table(metaDemog_uniq$meta.RecipientGender, metaDemog_uniq$meta.Outcome)
chisq.test(tbl)



tbl = table(metaDemog_uniq$meta.GenotypingClass, metaDemog_uniq$meta.Outcome)
tbl = tbl[-3,]
chisq.test(tbl)
