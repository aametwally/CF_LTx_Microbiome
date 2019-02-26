require("heatmap3")
library("data.table")
library(ggplot2)
library("phyloseq")
library("DESeq2")
require(vegan)
#library("Rtsne")
dev.off()
rm(list=ls())


setwd("/Users/ahmedmetwally/Dropbox/LungTransplantWASHU/CF_Paper_v4.0")
setwd("D:/Dropbox/LungTransplantWASHU/CF_Paper_v4.0")


# The palette with grey:
cbPalette  = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbbPalette = c("#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000", "#E69F00", "#56B4E9", "#009E73")
cbPalette  = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbbPalette_10 = c("#999999", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000", "#E69F00", "#56B4E9", "#009E73", "darkorchid")
cbbPalette_9 = c("#56B4E9", "sienna1", "#009E73", "darkorchid",  "#F0E442", "#0072B2", "#CC79A7", "#999999", "#000000")
cbbPalette_12 = c("#999999", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000", "#E69F00", "#56B4E9", "#009E73", "darkorchid",
                  "mistyrose1", "sienna1")





# Prepare metadata matrix
meta = read.csv(file="metadata/CF_12_v2.3.csv", header = TRUE)
View(meta)








###########################################
##### Correlate/Linear regression FEV1 with Phyla
###########################################
fit_bos <- lm(FEV1_Per_Predicted ~ Actinobacteria_Phyla_Percent, subset = Outcome == "BO",data = meta)
summary(fit_bos)
#p-value=0.0912, R2=0.11



fit_non_bos <- lm(FEV1_Per_Predicted ~ Actinobacteria_Phyla_Percent, subset = Outcome == "non-BO",data = meta)
summary(fit_non_bos)
#p-value=0.603, R2=0.006




## FACET
jpeg("LgTxCF12_FEV1_Actinobacteria_pergroup.jpg", res = 300, height = 8, width = 15, units = 'cm')
ggplot(meta, aes(x = Actinobacteria_Phyla_Percent, y = FEV1_Per_Predicted, color = Outcome)) + 
  geom_point() +
  stat_smooth(method = "lm", aes(color=Outcome)) +
  ggtitle("FEV (%) Perdicted vs Actinobacteria abundance") +  
  scale_fill_manual(values=c( "red", "midnightblue"),  breaks=c("BO", "non-BO")) + 
  scale_color_manual(values=c( "red", "midnightblue"),  breaks=c("BO", "non-BO")) +
  #scale_fill_manual(values = cbbPalette_12) +
  #scale_color_manual(values=cbbPalette_12) + 
  #scale_color_manual(values = cbbPalette_12) +
  #scale_fill_brewer(palette = "Set3")
  labs(x = "Actinobacteria abundance", y = "FEV (%) Perdicted") + 
  theme_bw() +
  theme(axis.text.x = element_text(colour="black", size=10, angle=0, hjust=0.5, vjust=0, face="bold"),
        axis.text.y = element_text(colour="black", size=10, angle=0, hjust=0.5, vjust=0.5, face="bold"),
        axis.title.x = element_text(colour="black", size=12, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=12, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=10, face="plain"), 
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
        legend.position="none",
        strip.text = element_text(size=12, face = "bold")) +
  guides(colour = guide_legend(ncol = 1)) + facet_grid(. ~ Outcome)
dev.off()




############## All without FACET
fit_all <- lm(FEV1_Per_Predicted ~ Actinobacteria_Phyla_Percent,data = meta)
summary(fit_all)

jpeg("LgTxCF12_FEV1_Actinobacteria_pergroup.jpg", res = 300, height = 8, width = 15, units = 'cm')
ggplot(meta, aes(x = Actinobacteria_Phyla_Percent, y = FEV1_Per_Predicted)) + 
  geom_point() +
  stat_smooth(method = "lm") +
  ggtitle("FEV (%) Perdicted vs Actinobacteria abundance") +  
  scale_fill_manual(values=c( "red", "midnightblue"),  breaks=c("BO", "non-BO")) + 
  scale_color_manual(values=c( "red", "midnightblue"),  breaks=c("BO", "non-BO")) +
  #scale_fill_manual(values = cbbPalette_12) +
  #scale_color_manual(values=cbbPalette_12) + 
  #scale_color_manual(values = cbbPalette_12) +
  #scale_fill_brewer(palette = "Set3")
  labs(x = "Actinobacteria abundance", y = "FEV (%) Perdicted") + 
  theme_bw() +
  theme(axis.text.x = element_text(colour="black", size=10, angle=0, hjust=0.5, vjust=0, face="bold"),
        axis.text.y = element_text(colour="black", size=10, angle=0, hjust=0.5, vjust=0.5, face="bold"),
        axis.title.x = element_text(colour="black", size=12, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=12, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=10, face="plain"), 
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
        legend.position="none",
        strip.text = element_text(size=12, face = "bold")) +
  guides(colour = guide_legend(ncol = 1))# + facet_grid(. ~ Outcome)
dev.off()





############################
############################
############################
meta$Proteobacteria_Phyla_Percent

fit_bos <- lm(FEV1_Per_Predicted ~ Proteobacteria_Phyla_Percent, subset = Outcome == "BO",data = meta)
summary(fit_bos)
#p-value=0.345, R2=0.03435



fit_non_bos <- lm(FEV1_Per_Predicted ~ Proteobacteria_Phyla_Percent, subset = Outcome == "non-BO",data = meta)
summary(fit_non_bos)
#p-value=0.0698, R2=0.07796

## FACET
jpeg("LgTxCF12_FEV1_Proteobacteria_pergroup.jpg", res = 300, height = 8, width = 15, units = 'cm')
ggplot(meta, aes(x = Proteobacteria_Phyla_Percent, y = FEV1_Per_Predicted, color = Outcome)) + 
  geom_point() +
  stat_smooth(method = "lm", aes(color=Outcome)) +
  ggtitle("FEV (%) Perdicted vs Proteobacteria abundance") +  
  scale_fill_manual(values=c( "red", "midnightblue"),  breaks=c("BO", "non-BO")) + 
  scale_color_manual(values=c( "red", "midnightblue"),  breaks=c("BO", "non-BO")) +
  #scale_fill_manual(values = cbbPalette_12) +
  #scale_color_manual(values=cbbPalette_12) + 
  #scale_color_manual(values = cbbPalette_12) +
  #scale_fill_brewer(palette = "Set3")
  labs(x = "Proteobacteria abundance", y = "FEV (%) Perdicted") + 
  theme_bw() +
  theme(axis.text.x = element_text(colour="black", size=10, angle=0, hjust=0.5, vjust=0, face="bold"),
        axis.text.y = element_text(colour="black", size=10, angle=0, hjust=0.5, vjust=0.5, face="bold"),
        axis.title.x = element_text(colour="black", size=12, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=12, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=10, face="plain"), 
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
        legend.position="none",
        strip.text = element_text(size=12, face = "bold")) +
  guides(colour = guide_legend(ncol = 1)) + facet_grid(. ~ Outcome)
dev.off()



###################################################
#### all bacteria  #######
###################################################
meta$Bacteroidetes_Phyla_Percent
meta$TimeAfterTransplantation

fit_bos <- lm(FEV1_Per_Predicted ~ Proteobacteria_Phyla_Percent + Actinobacteria_Phyla_Percent + Firmicutes_Phyla_Percent +
              Bacteroidetes_Phyla_Percent + TimeAfterTransplantation, subset = Outcome == "BO",data = meta)
summary(fit_bos)

# (Intercept)                  -34.66213   54.90404  -0.631   0.5341  
# Proteobacteria_Phyla_Percent  97.07444   50.87397   1.908   0.0689 .
# Actinobacteria_Phyla_Percent 750.49732  293.16783   2.560   0.0175 *
# Firmicutes_Phyla_Percent     254.69063  167.62732   1.519   0.1423  
# Bacteroidetes_Phyla_Percent         NA         NA      NA       NA  
# TimeAfterTransplantation      -0.01754    0.02679  -0.655   0.5192



fit_non_bos <- lm(FEV1_Per_Predicted ~ Proteobacteria_Phyla_Percent + Actinobacteria_Phyla_Percent + Firmicutes_Phyla_Percent +
                Bacteroidetes_Phyla_Percent + TimeAfterTransplantation, subset = Outcome == "non-BO",data = meta)
summary(fit_non_bos)
# Proteobacteria_Phyla_Percent  46.749479  50.532299   0.925    0.361
# Actinobacteria_Phyla_Percent  -4.323813 109.035409  -0.040    0.969
# Firmicutes_Phyla_Percent       4.882297 131.512766   0.037    0.971
# Bacteroidetes_Phyla_Percent          NA         NA      NA       NA
# TimeAfterTransplantation       0.008773   0.013227   0.663    0.511





###################################################
#### Subset samples based on MetaLonDA days #######
###################################################
meta_subset_200 = subset(meta, TimeAfterTransplantation>200)
dim(meta)
dim(meta_subset_200)
View(meta_subset_200)

fit_bos <- lm(FEV1_Per_Predicted ~ Proteobacteria_Phyla_Percent + Actinobacteria_Phyla_Percent + Firmicutes_Phyla_Percent +
                Bacteroidetes_Phyla_Percent + TimeAfterTransplantation, subset = Outcome == "BO",data = meta_subset_200)
summary(fit_bos)



fit_non_bos <- lm(FEV1_Per_Predicted ~ Firmicutes_Phyla_Percent + Bacteroidetes_Phyla_Percent + Proteobacteria_Phyla_Percent + Actinobacteria_Phyla_Percent, subset = Outcome == "non-BO",data = meta)
summary(fit_non_bos)





################################################
#### 
################################################