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

#library(MetaLonDA)
library(metagenomeSeq)
require("heatmap3")
library("data.table")
require(plyr)
#require(vegan)
library("phyloseq")
library("DESeq2")
# library("genefilter")
# library("ape")
library("phytools")
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



###############################
###############################
######### PhyloSeq ############
###############################
###############################
countTable = read.csv(file="taxProfiles/CF12/LungTxCF12_countMatrix.csv", header=TRUE, check.names = FALSE, row.names = 1)
taxaTable = read.csv(file="taxProfiles/CF12/LungTxCF12_OTUMatrix.csv", header=TRUE, row.names = 1)
meta = read.csv(file="metadata/CF_12_v1.3.csv", header = TRUE)
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




###############################################
##### Differential Abundance using DESeq #####
###############################################
### Function to calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
}

### Aglomorate
microbial_agglomartive = tax_glom(microbial_100, "Family")


dataset_deseq = microbial_agglomartive
microbial_dds = phyloseq_to_deseq2(dataset_deseq, ~ Outcome)
geoMeans = apply(counts(microbial_dds), 1, gm_mean)
microbial_dds_est = estimateSizeFactors(microbial_dds, geoMeans = geoMeans)
microbial_deseq = DESeq(microbial_dds_est, fitType="local")
res = results(microbial_deseq)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.05
sigtab = res[(res$pvalue < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(dataset_deseq)[rownames(sigtab), ], "matrix"))
head(sigtab)
write.csv(sigtab, file = "LgTxCF12_microbial_taxa_deseq_50.csv")
#posigtab = sigtab[sigtab[, "log2FoldChange"] > 0, ]
#posigtab = posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]

mcols(res, use.names = TRUE)
summary(res)

##############################
### PLot the deseq output ####
##############################
theme_set(theme_bw())
sigtabgen = subset(sigtab, !is.na(Genus))

# Phylum order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))

# Genus order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels=names(x))



jpeg("LgTxCF12_deseq_logfold.jpg", res = 300, height = 10, width = 20, units = 'cm')
ggplot(sigtabgen, aes(y=Genus, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=3) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
dev.off()





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

# source("/Users/ahmedmetwally/Dropbox/LungTransplantWASHU/v2.7_countries_cf12_package/R/Metalonda.R")
# source("/Users/ahmedmetwally/Dropbox/LungTransplantWASHU/v2.7_countries_cf12_package/R/CurveFitting.R")
# source("/Users/ahmedmetwally/Dropbox/LungTransplantWASHU/v2.7_countries_cf12_package/R/Permutation.R")
# source("/Users/ahmedmetwally/Dropbox/LungTransplantWASHU/v2.7_countries_cf12_package/R/Visualization.R")
# source("/Users/ahmedmetwally/Dropbox/LungTransplantWASHU/v2.7_countries_cf12_package/R/Normalization.R")


source("D:/Dropbox/metaLonDA_Package/v2.9_countries_cf12_package/R/Metalonda.R")
source("D:/Dropbox/metaLonDA_Package/v2.9_countries_cf12_package/R/CurveFitting.R")
source("D:/Dropbox/metaLonDA_Package/v2.9_countries_cf12_package/R/Permutation.R")
source("D:/Dropbox/metaLonDA_Package/v2.9_countries_cf12_package/R/Visualization.R")
source("D:/Dropbox/metaLonDA_Package/v2.9_countries_cf12_package/R/Normalization.R")






## TODO: reorder the metadata to follow the the order of the OTU table
## TODO: Glom one rank and replace the taxID with the actual name for more descriptive figures
data_real_count = as.data.frame(otu_table(microbial))
Group_real = as.vector(as.data.frame(sample_data(microbial))$Outcome)
ID_real = as.vector(as.data.frame(sample_data(microbial))$SUBJECT.WUTX)
Time_real = as.data.frame(sample_data(microbial))$TimeAfterTransplantation


output_all_nbinomial_2 = metalondaAll(Count = data_real_count[1:20,], Time = Time_real, Group = Group_real, 
                                         ID = ID_real, fit.method = "nbinomial", n.perm = 5, 
                                         num.intervals = 10, parall = FALSE, pvalue_threshold = 0.05, 
                                         adjust.method = "none", time.unit = "days", norm.method = "none",
                                      prefix = "Lx")

visualizeTimeIntervals(intervalDetails = output_all_nbinomial_2$output_summary, prefix = "LX")
write.csv(output_all_nbinomial_2$output_summary, file = "LgTxCF12_metalonda_css.csv", row.names = FALSE)



#### Test new metalonda
output_all_nbinomial_2 = metalondaAll(Count = data_real_count[1:20,], Time = Time_real, Group = Group_real, 
                                      ID = ID_real, fit.method = "nbinomial", n.perm = 5, 
                                      num.intervals = 99, parall = FALSE, pvalue.threshold = 0.05, 
                                      adjust.method = "none", time.unit = "days", norm.method = "none",
                                      prefix = "Lx")




### Test vector length
data_real_count
Group_real
ID_real
Time_real
Time_real = Time_real[-1]

if(length(Time_real) == length(ID_real))
{
  if(ncol(data_real_count) == length(Group_real))
  {
    if(length(Time_real) == length(Group_real))
    {
      cat("dimentionality check passed")
    } else
    {
      stop("The length of the annotation vectors don't match or not consistent 
           with the number of columns of the count matrix.")
    }
  } else
  {
    stop("The length of the annotation vectors don't match or not consistent 
           with the number of columns of the count matrix.")
  }
} else
{
  stop("The length of the annotation vectors don't match or not consistent 
           with the number of columns of the count matrix.")
}



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