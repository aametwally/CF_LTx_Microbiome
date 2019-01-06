#install.packages("MetaLonDA")
#install_github("aametwally/MetaLonDA", ref = "v1.1.2")

library("devtools")
library("MetaLonDA")
library("phyloseq")
library("DESeq2")
library("ggplot2")
library("zoo")

dev.off()
rm(list=ls())


setwd("/Users/ahmedmetwally/Dropbox/LungTransplantWASHU/CF_LTx_Microbiome/")
setwd("D:/Dropbox/LungTransplantWASHU/CF_LTx_Microbiome/")



# require(gss)
# require(plyr)
# require(caTools)
# require(ggplot2)
# require(doParallel)
# require(parallel)
# require(metagenomeSeq)
# require(DESeq2)
# require(edgeR)
# source("D:/Dropbox/metalona_dev/R/Metalonda.R")
# source("D:/Dropbox/metalona_dev/R/CurveFitting.R")
# source("D:/Dropbox/metalona_dev/R/Permutation.R")
# source("D:/Dropbox/metalona_dev/R/Visualization.R")
# source("D:/Dropbox/metalona_dev/R/Normalization.R")



###############################
###############################
######### PhyloSeq ############
###############################
###############################
countTable = read.csv(file="data/LungTxCF12_countMatrix.csv", header=TRUE, check.names = FALSE, row.names = 1)
taxaTable = read.csv(file="data/LungTxCF12_OTUMatrix.csv", header=TRUE, row.names = 1)
meta = read.csv(file="data/CF_12_v2.4.csv", header = TRUE)


## Prepare Phyloseq data object
OTU = otu_table(countTable, taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(taxaTable))
META = sample_data(meta)
sample_names(META) = meta$SampleIDFP
physeq = phyloseq(OTU, TAX, META)



##############################################
######## Normalization #######################
##############################################
## Relative abudnance
## microbiomeRA = transform_sample_counts(microbial, function(x) x / sum(x))

## Median of ratios:
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
}
phseq_dds = phyloseq_to_deseq2(physeq, ~ Outcome)
geoMeans = apply(counts(phseq_dds), 1, gm_mean)
phseq_dds_est = estimateSizeFactors(phseq_dds, geoMeans = geoMeans)
otu_matrix_norm = as.data.frame(counts(phseq_dds_est, normalized=TRUE))
OTU_norm = otu_table(otu_matrix_norm, taxa_are_rows = TRUE)
physeq_norm = phyloseq(OTU_norm, TAX, META)



#####################################################
############### Filteraion ##########################
#####################################################
## Filter out Eukaryota and root taxa
mic = subset_taxa(physeq_norm, Superkingdom != "Eukaryota")
microbial = subset_taxa(mic, Superkingdom != "")


### Retain taxa which appears in at leas 5% of samples with at minimum 5 reads. Based on: https://bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/MicrobiomeWorkflowII.html
level = "Species"
microbial.subset = subset_taxa(microbial, !is.na(level) & !level %in% c("", "uncharacterized"))
prevdf = apply(X = otu_table(microbial.subset),
               MARGIN = ifelse(taxa_are_rows(microbial.subset), yes = 1, no = 2),
               FUN = function(x){sum(x > 5)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(microbial.subset),
                    tax_table(microbial.subset))

plyr::ddply(prevdf, level, function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})


# Subset to the remaining phyla
prevdf1 = subset(prevdf, Species %in% get_taxa_unique(microbial.subset, "Species"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(microbial.subset),color=Species))+
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Species) + theme(legend.position="none")

prevalenceThreshold = 0.05 * nsamples(microbial.subset)
prevalenceThreshold
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
microbial.filtered = prune_taxa(keepTaxa, microbial.subset)




### Save Burkholderiaceae and Pseudomonadaceae families
Burkholderiaceae_Family = subset_taxa(microbial.filtered, Family == "Burkholderiaceae")
save(Burkholderiaceae_Family, file = "Burkholderiaceae_Family_PhyloseqObject.RData")

Pseudomonadaceae_Family = subset_taxa(microbial.filtered, Family == "Pseudomonadaceae")
save(Pseudomonadaceae_Family, file = "Pseudomonadaceae_Family_PhyloseqObject.RData")

View(otu_table(Pseudomonadaceae_Family))
View(tax_table(Pseudomonadaceae_Family))




#####################################################
############### MetaLonDA ###########################
#####################################################
## Change level to change the rank on which MetaLonDA needs to test 
level =  "Phylum"
microbial.glom = tax_glom(microbial.filtered, level)
dim(otu_table(microbial.glom))

apply(otu_table(microbial.glom), 1 , sum)
apply(otu_table(microbial.glom), 1 , mean)
apply(otu_table(microbial.glom), 1 , median)

## Rename rownames to be bacteria microbial name instead of the TaxID
microbial.glom.count = as.data.frame(otu_table(microbial.glom))
for (i in rownames(microbial.glom.count))
{
  x = as.vector(tax_table(microbial)[which(rownames(tax_table(microbial)) == i), level])
  lab = sprintf("%s_%s", x, i)
  rownames(microbial.glom.count)[which(rownames(microbial.glom.count) == i)] = lab
  cat(i, "\n")
}


### Execlude from the analysis
# which(rownames(tmp) == "Thermodesulfobacteria_1295609")
# tmp = tmp[-18, ]

Group_real = as.vector(as.data.frame(sample_data(microbial.glom))$Outcome)
ID_real = as.vector(as.data.frame(sample_data(microbial.glom))$AnnotatedID)
Time_real = as.data.frame(sample_data(microbial.glom))$TimeAfterTransplantation

output_all_nbinomial_2 = metalondaAll(Count = microbial.glom.count, Time = Time_real, Group = Group_real, 
                                      ID = ID_real, fit.method = "nbinomial", n.perm = 1000, 
                                      num.intervals = 99, parall = FALSE, pvalue.threshold = 0.1, 
                                      adjust.method = "BH", time.unit = "days", norm.method = "none",
                                      prefix = "LTx_Phylum", col = c("firebrick", "blue"))





########################################
########## MetaLonDA for Burlkerderia (all ranks) ##
########################################

Burk_Order = subset_taxa(microbial.filtered, Order == "Burkholderiales")
############### Agglomerate Taxa ####################

level =  "Order"
Burk_Order.glom = tax_glom(Burk_Order, level)
dim(otu_table(Burk_Order.glom))


level =  "Family"
Burk_Order.glom_family = tax_glom(Burk_Order, level)
dim(otu_table(Burk_Order.glom_family))


level =  "Genus"
Burk_Order.glom_genus = tax_glom(Burk_Order, level)
dim(otu_table(Burk_Order.glom_genus))

level =  "Species"
Burk_Order.glom_spices = tax_glom(Burk_Order, level)
dim(otu_table(Burk_Order.glom_spices))


### MetaLonDA for all burckolderia (genus) species
## Rename rownames to be bacteria microbial name instead of the TaxID
## TODO: replace "Burk_Order.glom_spices" argument with the rank that wanted to be tested
microbial.glom = Burk_Order.glom_spices ## TODO: Change this if you want to change the rank
level = "Species"  ### TODO: Change this if you want to change the rank
microbial.glom.count = as.data.frame(otu_table(microbial.glom))
for (i in rownames(microbial.glom.count))
{
  x = as.vector(tax_table(microbial)[which(rownames(tax_table(microbial)) == i), level])
  lab = sprintf("%s_%s", x, i)
  rownames(microbial.glom.count)[which(rownames(microbial.glom.count) == i)] = lab
  cat(i, "\n")
}


### Exclude from the analysis
# which(rownames(tmp) == "Thermodesulfobacteria_1295609")
# tmp = tmp[-18, ]

Group_real = as.vector(as.data.frame(sample_data(microbial.glom))$Outcome)
ID_real = as.vector(as.data.frame(sample_data(microbial.glom))$AnnotatedID)
Time_real = as.data.frame(sample_data(microbial.glom))$TimeAfterTransplantation

output_all_nbinomial_2 = metalondaAll(Count = microbial.glom.count, Time = Time_real, Group = Group_real, 
                                      ID = ID_real, fit.method = "nbinomial", n.perm = 1000, 
                                      num.intervals = 99, parall = FALSE, pvalue.threshold = 0.1, 
                                      adjust.method = "BH", time.unit = "days", norm.method = "none",
                                      prefix = "Lx", col = c("firebrick", "blue"))











#####################################################
############### OmicsLonDA ###########################
#####################################################
library(ggplot2)
library(gss)
library(plyr)
library(pracma)
library(parallel)
library(doParallel)
library(zoo)
library(pheatmap)
library(viridis)
source("/Users/ahmedmetwally/Dropbox/metalonda_work/OmicsLonDA_dev/OmicsLonDA/R/OmicsLonDA.R")
source("/Users/ahmedmetwally/Dropbox/metalonda_work/OmicsLonDA_dev/OmicsLonDA/R/CurveFitting.R")
source("/Users/ahmedmetwally/Dropbox/metalonda_work/OmicsLonDA_dev/OmicsLonDA/R/Visualization.R")
source("/Users/ahmedmetwally/Dropbox/metalonda_work/OmicsLonDA_dev/OmicsLonDA/R/Permutation.R")
source("/Users/ahmedmetwally/Dropbox/metalonda_work/OmicsLonDA_dev/OmicsLonDA/R/Normalization.R")


count.table = microbial.glom.count
metadata.table = sample_data(microbial.glom)
colnames(metadata.table)[which(colnames(metadata.table) == "Outcome")] = "Group"
colnames(metadata.table)[which(colnames(metadata.table) == "AnnotatedID")] = "Subject"
colnames(metadata.table)[which(colnames(metadata.table) == "TimeAfterTransplantation")] = "Time"
metadata.table = metadata.table[,c("Time", "Group", "Subject")]
metadata.table$Group = as.vector(metadata.table$Group)
metadata.table$Subjet = as.vector(metadata.table$Subject)


output.omicslondaAll = omicslondaAll(formula = Count ~ Time, countTable = count.table, metadata = metadata.table, 
                                     n.perm = 1000, fit.method = "ssgaussian", num.intervals = 99,
                                     parall = FALSE, pvalue.threshold = 0.1,     
                                     adjust.method = "BH", time.unit = "hours", norm.method = "none", prefix = "OmicsLonDA_Phylum", 
                                     ylabel = "Normalized Read Count", col = c("firebrick", "blue"))
