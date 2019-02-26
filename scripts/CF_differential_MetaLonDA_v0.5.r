#install.packages("MetaLonDA")
library(devtools)
install_github("aametwally/MetaLonDA", ref = "v1.1.2")
library(MetaLonDA)
library("phyloseq")
library(DESeq2)
library("ggplot2")
dev.off()
rm(list=ls())
setwd("/Users/ahmedmetwally/Dropbox/LungTransplantWASHU/CF_Paper_v4.0")
setwd("D:/Dropbox/LungTransplantWASHU/CF_Paper_v4.0")



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
countTable = read.csv(file="../taxProfiles/CF12/LungTxCF12_countMatrix.csv", header=TRUE, check.names = FALSE, row.names = 1)
taxaTable = read.csv(file="../taxProfiles/CF12/LungTxCF12_OTUMatrix.csv", header=TRUE, row.names = 1)
meta = read.csv(file="metadata/CF_12_v2.3.csv", header = TRUE)
# meta$SUBJECT.WUTX = as.factor(meta$SUBJECT.WUTX)
# meta$MiSeqRunNum = as.factor(meta$MiSeqRunNum)
# class(meta$Outcome)


## Prepare Phyloseq data object
OTU = otu_table(countTable, taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(taxaTable))
META = sample_data(meta)
sample_names(META) = meta$SampleIDFP
physeq = phyloseq(OTU, TAX, META)



##############################################
######## Normalization #######################
##############################################
## TODO: Try different normalization methods

## RA
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

View(otu_table(physeq_norm))
View(otu_table(physeq))
#tmp = otu_matrix_norm

#####################################################
############### Filteraion ##########################
#####################################################
# microbiomeRA = transform_sample_counts(microbial, function(x) x / sum(x))
#bacteria = subset_taxa(microbiomeRA, Superkingdom == "Bacteria")
#virus = subset_taxa(microbiomeRA, Superkingdom == "Viruses")
#archaea = subset_taxa(microbiomeRA, Superkingdom = "Archaea")
## Filter based on variance
# hist(log10(apply(otu_table(microbial), 1, var)), xlab = "log10(variance)", main = "A large fraction of OTUs have very low variance")
# varianceThreshold = 30
# keepOTUs = apply(otu_table(microbial), 1, var) > varianceThreshold
# microbial_varfiltered = prune_taxa(keepOTUs, microbial)
#microbial = subset_taxa(microbial, Superkingdom != "Viruses")
#unique(sort(tax_table(microbial)[,"Superkingdom"]))




## Filter out Eukaryota and root taxa
mic = subset_taxa(physeq_norm, Superkingdom != "Eukaryota")
microbial = subset_taxa(mic, Superkingdom != "")


### Tutorial: https://bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/MicrobiomeWorkflowII.html
#level =  "Phylum"
#microbial.glom = tax_glom(microbial, level)




### Retain taxa which appears in at leas 5% of samples with at minimum 5 reads.
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



microbial.filtered
microbial



#####################################################
############### Agglomerate Taxa ####################
#####################################################
microbial.filtered
level =  "Family"
microbial.glom = tax_glom(microbial.filtered, level)
dim(otu_table(microbial.glom))
#View(tax_table(microbial.glom))
#View(otu_table(microbial.glom))


apply(otu_table(microbial.glom), 1 , sum)
apply(otu_table(microbial.glom), 1 , mean)
apply(otu_table(microbial.glom), 1 , median)



#####################################################
############### MetaLonDA ###########################
#####################################################
microbial.glom

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
                                      num.intervals = 99, parall = FALSE, pvalue.threshold = 0.05, 
                                      adjust.method = "none", time.unit = "days", norm.method = "none",
                                      prefix = "Lx", col = c("firebrick", "blue"))




