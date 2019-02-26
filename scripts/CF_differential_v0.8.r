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

# install.packages("MetaLonDA")
# library(MetaLonDA)
library(metagenomeSeq)
require("heatmap3")
library("data.table")
require(plyr)
#require(vegan)
library("phyloseq")
library("DESeq2")
# library("genefilter")
# library("ape")
#library("phytools")
# library("vegan")
library("ggplot2")
# library("plot3D")
# library("calibrate")
# library("ggbiplot")
# library(reshape)
# library(ggplot2)
# library(limma)
#require(cowplot)
# library(gridGraphics)
dev.off()
rm(list=ls())


setwd("/Users/ahmedmetwally/Dropbox/LungTransplantWASHU")
setwd("D:/Dropbox/LungTransplantWASHU")


#############
### TODO: row.names = 1 in all other files
#####


##TODO: Reinstall DESEQ2
## USe the iterative method for geomtric mean
# use lfc shrick before drawing the MA plot
# packageVersion("DESeq2")

###############################
###############################
######### PhyloSeq ############
###############################
###############################
countTable = read.csv(file="taxProfiles/CF12/LungTxCF12_countMatrix.csv", header=TRUE, check.names = FALSE, row.names = 1)
taxaTable = read.csv(file="taxProfiles/CF12/LungTxCF12_OTUMatrix.csv", header=TRUE, row.names = 1)
meta = read.csv(file="metadata/CF_12_v1.4.csv", header = TRUE)
# meta$SUBJECT.WUTX = as.factor(meta$SUBJECT.WUTX)
# meta$MiSeqRunNum = as.factor(meta$MiSeqRunNum)
# class(meta$Outcome)

# #### Assign samples to bins
# meta$bin = 0
# meta$bin[meta$TimeAfterTransplantation<50] = 50
# meta$bin[abs(meta$TimeAfterTransplantation - 100) < 50] = 100
# meta$bin[abs(meta$TimeAfterTransplantation - 200) < 50] = 200
# meta$bin[abs(meta$TimeAfterTransplantation - 300) < 50] = 300
# meta$bin[abs(meta$TimeAfterTransplantation - 400) < 50] = 400
# meta$bin[abs(meta$TimeAfterTransplantation - 500) < 50 | meta$TimeAfterTransplantation>500] = 500


## Prepare Phyloseq data object
OTU = otu_table(countTable, taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(taxaTable))
META = sample_data(meta)
sample_names(META) = meta$SampleIDFP
physeq = phyloseq(OTU, TAX, META)


## Filter out Eukaryota and root taxa
mic = subset_taxa(physeq, Superkingdom != "Eukaryota")
microbial = subset_taxa(mic, Superkingdom != "")
microbial = subset_taxa(microbial, Superkingdom != "Viruses")
#unique(sort(tax_table(microbial)[,"Superkingdom"]))


##############################################
######## Normalization #######################
##############################################
## TODO: Try different normalization methods
# microbiomeRA = transform_sample_counts(microbial, function(x) x / sum(x))

#bacteria = subset_taxa(microbiomeRA, Superkingdom == "Bacteria")
#virus = subset_taxa(microbiomeRA, Superkingdom == "Viruses")
#archaea = subset_taxa(microbiomeRA, Superkingdom = "Archaea")



##################################
#### Cross sectional analysis ####
##################################
microbial_50 <- subset_samples(microbial, bin==50)
microbial_100 <- subset_samples(microbial, bin==100)
microbial_200 <- subset_samples(microbial, bin==200)
microbial_300 <- subset_samples(microbial, bin==300)
microbial_400 <- subset_samples(microbial, bin==400)
microbial_500 <- subset_samples(microbial, bin==500)
microbial_last <- subset_samples(microbial, TimePoint_Status=="last")




## TODO: Take out all viruses
## TODO: change the pvalue to adj
## TODO: Take out the low variance species
hist(log10(apply(otu_table(microbial), 1, var)), xlab = "log10(variance)", main = "A large fraction of OTUs have very low variance")
varianceThreshold = 30
keepOTUs = apply(otu_table(microbial), 1, var) > varianceThreshold
microbial_varfiltered = prune_taxa(keepOTUs, microbial)

###############################################
##### Differential Abundance using DESeq #####
###############################################
### Function to calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
}

### Aglomorate
## TODO: Check if the total is based on one level or the aglom is wrong
## TODO: Check which level is over abundant or lower abundant
# results will extract the results table for a comparison of the last level over the first level. 
# res <- results(dds, contrast=c("dex","trt","untrt"))
# mcols(res, use.names = TRUE)
# non-bos vs bos
# # summary(res)

######################################
### Do it for all samples (Outcome) ##
######################################
dataset_deseq = microbial
microbial_dds = phyloseq_to_deseq2(dataset_deseq, ~ Outcome)
geoMeans = apply(counts(microbial_dds), 1, gm_mean)
microbial_dds_est = estimateSizeFactors(microbial_dds, geoMeans = geoMeans)
microbial_deseq = DESeq(microbial_dds_est, fitType="local")
res = results(microbial_deseq, contrast=c("Outcome","bos","non-bos"))
res = res[order(res$padj, na.last=NA), ]
alpha = 0.05
sigtab = res[(res$pvalue < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(dataset_deseq)[rownames(sigtab), ], "matrix"))
head(sigtab)
View(sigtab)
write.csv(res, file = "LgTxCF12_microbial_taxa_deseq_microbial_all.csv")




######################################
### Do it for all samples (Outcome): Phylum ##
######################################
microbial_agglomartive = tax_glom(microbial, "Phylum")
dataset_deseq = microbial_agglomartive
microbial_dds = phyloseq_to_deseq2(dataset_deseq, ~ Outcome)
geoMeans = apply(counts(microbial_dds), 1, gm_mean)
microbial_dds_est = estimateSizeFactors(microbial_dds, geoMeans = geoMeans)
microbial_deseq = DESeq(microbial_dds_est, fitType="local")
res = results(microbial_deseq, contrast=c("Outcome","bos","non-bos"))
res = res[order(res$padj, na.last=NA), ]
alpha = 0.05
sigtab = res[(res$pvalue < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(dataset_deseq)[rownames(sigtab), ], "matrix"))
head(sigtab)
View(sigtab)
write.csv(res, file = "LgTxCF12_microbial_taxa_deseq_microbial_all_phylum.csv")



######################################
### Do it for all samples (Outcome): Family ##
######################################
microbial_agglomartive = tax_glom(microbial, "Family")
dataset_deseq = microbial_agglomartive
microbial_dds = phyloseq_to_deseq2(dataset_deseq, ~ Outcome)
geoMeans = apply(counts(microbial_dds), 1, gm_mean)
microbial_dds_est = estimateSizeFactors(microbial_dds, geoMeans = geoMeans)
microbial_deseq = DESeq(microbial_dds_est, fitType="local")
res = results(microbial_deseq, contrast=c("Outcome","bos","non-bos"))
res = res[order(res$padj, na.last=NA), ]
alpha = 0.05
sigtab = res[(res$pvalue < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(dataset_deseq)[rownames(sigtab), ], "matrix"))
head(sigtab)
View(sigtab)
write.csv(res, file = "LgTxCF12_microbial_taxa_deseq_microbial_all_family.csv")



######################################
### Do it for all samples (Outcome): Genus ##
######################################
microbial_agglomartive = tax_glom(microbial, "Genus")
dataset_deseq = microbial_agglomartive
microbial_dds = phyloseq_to_deseq2(dataset_deseq, ~ Outcome)
geoMeans = apply(counts(microbial_dds), 1, gm_mean)
microbial_dds_est = estimateSizeFactors(microbial_dds, geoMeans = geoMeans)
microbial_deseq = DESeq(microbial_dds_est, fitType="local")
res = results(microbial_deseq, contrast=c("Outcome","bos","non-bos"))
res = res[order(res$padj, na.last=NA), ]
alpha = 0.05
sigtab = res[(res$pvalue < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(dataset_deseq)[rownames(sigtab), ], "matrix"))
head(sigtab)
View(sigtab)
write.csv(res, file = "LgTxCF12_microbial_taxa_deseq_microbial_all_genus.csv")


######################################
### Do it for all samples (Outcome): Species ##
######################################
microbial_agglomartive = tax_glom(microbial, "Species")
dataset_deseq = microbial_agglomartive
microbial_dds = phyloseq_to_deseq2(dataset_deseq, ~ Outcome)
geoMeans = apply(counts(microbial_dds), 1, gm_mean)
microbial_dds_est = estimateSizeFactors(microbial_dds, geoMeans = geoMeans)
microbial_deseq = DESeq(microbial_dds_est, fitType="local")
res = results(microbial_deseq, contrast=c("Outcome","bos","non-bos"))
res = res[order(res$padj, na.last=NA), ]
alpha = 0.05
sigtab = res[(res$pvalue < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(dataset_deseq)[rownames(sigtab), ], "matrix"))
head(sigtab)
View(sigtab)
write.csv(res, file = "LgTxCF12_microbial_taxa_deseq_microbial_all_species.csv")


#######################
######  LAST ##########
#######################

######################################
### Do it for all samples: LAST (Outcome) ##
######################################
dataset_deseq = microbial_last
microbial_dds = phyloseq_to_deseq2(dataset_deseq, ~ Outcome)
geoMeans = apply(counts(microbial_dds), 1, gm_mean)
microbial_dds_est = estimateSizeFactors(microbial_dds, geoMeans = geoMeans)
microbial_deseq = DESeq(microbial_dds_est, fitType="local")
res = results(microbial_deseq, contrast=c("Outcome","bos","non-bos"))
res = res[order(res$padj, na.last=NA), ]
alpha = 0.05
sigtab = res[(res$pvalue < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(dataset_deseq)[rownames(sigtab), ], "matrix"))
head(sigtab)
View(sigtab)
write.csv(res, file = "LgTxCF12_microbial_taxa_deseq_microbial_last.csv")
sigtab_last_all = sigtab



######################################
### Do it for all samples (Outcome): : LAST: Phylum ##
######################################
microbial_agglomartive = tax_glom(microbial_last, "Phylum")
dataset_deseq = microbial_agglomartive
microbial_dds = phyloseq_to_deseq2(dataset_deseq, ~ Outcome)
geoMeans = apply(counts(microbial_dds), 1, gm_mean)
microbial_dds_est = estimateSizeFactors(microbial_dds, geoMeans = geoMeans)
microbial_deseq = DESeq(microbial_dds_est, fitType="local")
res = results(microbial_deseq, contrast=c("Outcome","bos","non-bos"))
res = res[order(res$padj, na.last=NA), ]
alpha = 0.05
sigtab = res[(res$pvalue < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(dataset_deseq)[rownames(sigtab), ], "matrix"))
head(sigtab)
View(sigtab)
write.csv(res, file = "LgTxCF12_microbial_taxa_deseq_microbial_last_phylum.csv")

sigtab_last_phylum = sigtab

######################################
### Do it for all samples (Outcome): LAST: Family ##
######################################
microbial_agglomartive = tax_glom(microbial_last, "Family")
dataset_deseq = microbial_agglomartive
microbial_dds = phyloseq_to_deseq2(dataset_deseq, ~ Outcome)
geoMeans = apply(counts(microbial_dds), 1, gm_mean)
microbial_dds_est = estimateSizeFactors(microbial_dds, geoMeans = geoMeans)
microbial_deseq = DESeq(microbial_dds_est, fitType="local")
res = results(microbial_deseq, contrast=c("Outcome","bos","non-bos"))
res = res[order(res$padj, na.last=NA), ]
alpha = 0.05
sigtab = res[(res$pvalue < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(dataset_deseq)[rownames(sigtab), ], "matrix"))
head(sigtab)
View(sigtab)
write.csv(res, file = "LgTxCF12_microbial_taxa_deseq_microbial_last_family.csv")
sigtab_last_family = sigtab


######################################
### Do it for all samples (Outcome): LAST: Genus ##
######################################
microbial_agglomartive = tax_glom(microbial_last, "Genus")
dataset_deseq = microbial_agglomartive
microbial_dds = phyloseq_to_deseq2(dataset_deseq, ~ Outcome)
geoMeans = apply(counts(microbial_dds), 1, gm_mean)
microbial_dds_est = estimateSizeFactors(microbial_dds, geoMeans = geoMeans)
microbial_deseq = DESeq(microbial_dds_est, fitType="local")
res = results(microbial_deseq, contrast=c("Outcome","bos","non-bos"))
res = res[order(res$padj, na.last=NA), ]
alpha = 0.05
sigtab = res[(res$pvalue < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(dataset_deseq)[rownames(sigtab), ], "matrix"))
head(sigtab)
View(sigtab)
write.csv(res, file = "LgTxCF12_microbial_taxa_deseq_microbial_last_genus.csv")
sigtab_last_genus = sigtab

######################################
### Do it for all samples (Outcome): LAST Species ##
######################################
microbial_agglomartive = tax_glom(microbial_last, "Species")
dataset_deseq = microbial_agglomartive
microbial_dds = phyloseq_to_deseq2(dataset_deseq, ~ Outcome)
geoMeans = apply(counts(microbial_dds), 1, gm_mean)
microbial_dds_est = estimateSizeFactors(microbial_dds, geoMeans = geoMeans)
microbial_deseq = DESeq(microbial_dds_est, fitType="local")
res = results(microbial_deseq, contrast=c("Outcome","bos","non-bos"))
res = res[order(res$padj, na.last=NA), ]
alpha = 0.05
sigtab = res[(res$pvalue < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(dataset_deseq)[rownames(sigtab), ], "matrix"))
head(sigtab)
View(sigtab)
write.csv(res, file = "LgTxCF12_microbial_taxa_deseq_microbial_last_species.csv")
sigtab_last_species = sigtab










#####################
###### 50  ##########
#####################


######################################
### Do it for all samples: 50 (Outcome) ##
######################################
dataset_deseq = microbial_50
microbial_dds = phyloseq_to_deseq2(dataset_deseq, ~ Outcome)
geoMeans = apply(counts(microbial_dds), 1, gm_mean)
microbial_dds_est = estimateSizeFactors(microbial_dds, geoMeans = geoMeans)
microbial_deseq = DESeq(microbial_dds_est, fitType="local")
res = results(microbial_deseq, contrast=c("Outcome","bos","non-bos"))
res = res[order(res$padj, na.last=NA), ]
alpha = 0.05
sigtab = res[(res$pvalue < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(dataset_deseq)[rownames(sigtab), ], "matrix"))
head(sigtab)
View(sigtab)
write.csv(res, file = "LgTxCF12_microbial_taxa_deseq_microbial_50.csv")




######################################
### Do it for all samples (Outcome): : 50: Phylum ##
######################################
microbial_agglomartive = tax_glom(microbial_50, "Phylum")
dataset_deseq = microbial_agglomartive
microbial_dds = phyloseq_to_deseq2(dataset_deseq, ~ Outcome)
geoMeans = apply(counts(microbial_dds), 1, gm_mean)
microbial_dds_est = estimateSizeFactors(microbial_dds, geoMeans = geoMeans)
microbial_deseq = DESeq(microbial_dds_est, fitType="local")
res = results(microbial_deseq, contrast=c("Outcome","bos","non-bos"))
res = res[order(res$padj, na.last=NA), ]
alpha = 0.05
sigtab = res[(res$pvalue < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(dataset_deseq)[rownames(sigtab), ], "matrix"))
head(sigtab)
View(sigtab)
write.csv(res, file = "LgTxCF12_microbial_taxa_deseq_microbial_50_phylum.csv")



######################################
### Do it for all samples (Outcome): 50 Family ##
######################################
microbial_agglomartive = tax_glom(microbial_50, "Family")
dataset_deseq = microbial_agglomartive
microbial_dds = phyloseq_to_deseq2(dataset_deseq, ~ Outcome)
geoMeans = apply(counts(microbial_dds), 1, gm_mean)
microbial_dds_est = estimateSizeFactors(microbial_dds, geoMeans = geoMeans)
microbial_deseq = DESeq(microbial_dds_est, fitType="local")
res = results(microbial_deseq, contrast=c("Outcome","bos","non-bos"))
res = res[order(res$padj, na.last=NA), ]
alpha = 0.05
sigtab = res[(res$pvalue < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(dataset_deseq)[rownames(sigtab), ], "matrix"))
head(sigtab)
View(sigtab)
write.csv(res, file = "LgTxCF12_microbial_taxa_deseq_microbial_50_family.csv")



######################################
### Do it for all samples (Outcome): 50 Genus ##
######################################
microbial_agglomartive = tax_glom(microbial_50, "Genus")
dataset_deseq = microbial_agglomartive
microbial_dds = phyloseq_to_deseq2(dataset_deseq, ~ Outcome)
geoMeans = apply(counts(microbial_dds), 1, gm_mean)
microbial_dds_est = estimateSizeFactors(microbial_dds, geoMeans = geoMeans)
microbial_deseq = DESeq(microbial_dds_est, fitType="local")
res = results(microbial_deseq, contrast=c("Outcome","bos","non-bos"))
res = res[order(res$padj, na.last=NA), ]
alpha = 0.05
sigtab = res[(res$pvalue < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(dataset_deseq)[rownames(sigtab), ], "matrix"))
head(sigtab)
View(sigtab)
write.csv(res, file = "LgTxCF12_microbial_taxa_deseq_microbial_50_genus.csv")


######################################
### Do it for all samples (Outcome): 50: Species ##
######################################
microbial_agglomartive = tax_glom(microbial_50, "Species")
dataset_deseq = microbial_agglomartive
microbial_dds = phyloseq_to_deseq2(dataset_deseq, ~ Outcome)
geoMeans = apply(counts(microbial_dds), 1, gm_mean)
microbial_dds_est = estimateSizeFactors(microbial_dds, geoMeans = geoMeans)
microbial_deseq = DESeq(microbial_dds_est, fitType="local")
res = results(microbial_deseq, contrast=c("Outcome","bos","non-bos"))
res = res[order(res$padj, na.last=NA), ]
alpha = 0.05
sigtab = res[(res$pvalue < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(dataset_deseq)[rownames(sigtab), ], "matrix"))
head(sigtab)
View(sigtab)
write.csv(res, file = "LgTxCF12_microbial_taxa_deseq_microbial_50_species.csv")






##################################################
###### Visualize differentially abundant taxa  ###
##################################################
# Species order
sigtabgen = subset(sigtab, !is.na(Species))
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Species, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Species = factor(as.character(sigtabgen$Species), levels=names(x))

sigtabgen$enriched[which(sigtabgen$log2FoldChange>0)] = "bos"
sigtabgen$enriched[which(sigtabgen$log2FoldChange<0)] = "non-bos"

## gather 
sigtabgen_selected = sigtabgen[,c("log2FoldChange", "pvalue", "padj", "Species", "enriched")]
colnames(sigtabgen_selected)[4] = "rank"


jpeg("LgTxCF12_differentiallyAbundant_taxa.jpg", res = 300, height = 15, width = 25, units = 'cm')
ggplot(sigtabgen_selected, aes(y=rank, x=log2FoldChange, color=enriched)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=3) +  
  ggtitle("Differentially abundant taxa (p-adj<0.05)") +  
  scale_fill_manual(values=c( "red", "midnightblue"),  breaks=c("bos", "non-bos")) + 
  scale_color_manual(values=c( "red", "midnightblue"),  breaks=c("bos", "non-bos")) +
  #scale_fill_brewer(palette = "Set3")
  labs( y = "Taxa", x = "Log2FoldChange") + 
  theme_bw() + #stat_summary(fun.y=mean, geom="point", size=1, color="black") + geom_boxplot(width=0.1) +
  #scale_y_continuous(breaks = c(0, 1)) + 
  theme(axis.text.x = element_text(colour="black", size=10, angle=0, hjust=0.5, vjust=0.5, face="bold"),
        axis.text.y = element_text(colour="black", size=10, angle=0, hjust=1, vjust=0.5, face="bold"),
        axis.title.x = element_text(colour="black", size=12, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=12, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=15, face="plain"), 
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position="top",
        strip.text = element_text(size=12, face = "bold")) +
  guides(colour = guide_legend(nrow = 1))
dev.off()







## Integrate all last bins

sigtab_last_all
sigtab_last_phylum
sigtab_last_family
sigtab_last_genus
sigtab_last_species










############################
##### Test Sparsity:  ######
############################
# TODO: “On the other hand, very sparse count datasets, with large counts for single samples per row and the 
# rest at 0, don’t fit well to the negative binomial distribution. Here, the VST or simply shifted log,  
# log(count+k), might be a safer choice than the rlog. A way that I test for sparsity is looking at 
# a plot of the row sum of counts and the proportion of count which is in a single sample.”

reads_sum <- rowSums(counts(microbial_dds_est))
reads_max <- apply(counts(microbial_dds_est), 1, max)
plot(reads_sum+1, reads_max/reads_sum, log="x")



##############################
### PLot the deseq output ####
##############################
### MA plots
res <- lfcShrink(microbial_dds_est,  res=res)
plotMA(res)



### TODO: Plot the one above with facet and for all variables
































#####################################################
############### MetaLonDA ###########################
#####################################################
#library(MetaLonDA)
require(gss)
require(plyr)
require(caTools)
require(ggplot2)
require(doParallel)
require(parallel)
require(metagenomeSeq)
require(DESeq2)
require(edgeR)


source("D:/Dropbox/metalona_dev/R/Metalonda.R")
source("D:/Dropbox/metalona_dev/R/CurveFitting.R")
source("D:/Dropbox/metalona_dev/R/Permutation.R")
source("D:/Dropbox/metalona_dev/R/Visualization.R")
source("D:/Dropbox/metalona_dev/R/Normalization.R")



#### Test metalonda




### Glom
microbial
level =  "Genus"
microbial.glom = tax_glom(microbial, level)
dim(otu_table(microbial.glom))
View(tax_table(microbial.glom))

## Normalization
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
}
phseq_dds = phyloseq_to_deseq2(microbial.glom, ~ Outcome)
geoMeans = apply(counts(phseq_dds), 1, gm_mean)
phseq_dds_est = estimateSizeFactors(phseq_dds, geoMeans = geoMeans)
otu_matrix_norm = as.data.frame(counts(phseq_dds_est, normalized=TRUE))
tmp = otu_matrix_norm


## Rename rownames to be bacteria microbial name instead of the TaxID
for (i in rownames(tmp))
{
  x = as.vector(tax_table(microbial)[which(rownames(tax_table(microbial)) == i), level])
  lab = sprintf("%s_%s", x, i)
  rownames(tmp)[which(rownames(tmp) == i)] = lab #as.vector(tax_table(microbial)[which(rownames(tax_table(microbial)) == i), "Species"])
  cat(i, "\n")
}

### Execlude from the analysis
# which(rownames(tmp) == "Thermodesulfobacteria_1295609")
# tmp = tmp[-18, ]


Group_real = as.vector(as.data.frame(sample_data(microbial))$Outcome)
ID_real = as.vector(as.data.frame(sample_data(microbial))$SUBJECT.WUTX)
Time_real = as.data.frame(sample_data(microbial))$TimeAfterTransplantation

output_all_nbinomial_2 = metalondaAll(Count = tmp, Time = Time_real, Group = Group_real, 
                                      ID = ID_real, fit.method = "nbinomial", n.perm = 100, 
                                      num.intervals = 99, parall = FALSE, pvalue.threshold = 0.05, 
                                      adjust.method = "BH", time.unit = "days", norm.method = "none",
                                      prefix = "Lx")


###############################
##### GLM ##########
###############################
tx = glm(Outcome ~ FEV1_beforeTx, data=meta, family = binomial)
summary(tx)
anova(tx)


tx_time = glm(Outcome ~ FVC_predicted, data = meta, family = binomial)
summary(tx_time)
anova(tx_time)
1-pchisq(96.61,72)
plot(tx_time)




##### Draw the time-series of Psuedomonas aeruginosa 
aeruginosa