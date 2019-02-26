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

install.packages("MetaLonDA")
library(MetaLonDA)
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



###########################
###### Network  ###########
###########################









