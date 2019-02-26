# install.packages("ggbiplot")
# install.packages("heatmap3")
# install.packages("phyloseq")
# install.packages("vegan")
# install.packages("DESeq2")
# install.packages("phyloseq")
# #install.packages("Rtsne")

# source("https://bioconductor.org/biocLite.R")
# biocLite("phyloseq")
# biocLite("ade4")

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
setwd("D:/Dropbox/LungTransplantWASHU/CF_Paper_v4.1")


# The palette with grey:
cbPalette  = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbbPalette = c("#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000", "#E69F00", "#56B4E9", "#009E73")
cbPalette  = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbbPalette_10 = c("#999999", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000", "#E69F00", "#56B4E9", "#009E73", "darkorchid")
cbbPalette_9 = c("#56B4E9", "sienna1", "#009E73", "darkorchid",  "#F0E442", "#0072B2", "#CC79A7", "#999999", "#000000")
cbbPalette_12 = c("#999999", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000", "#E69F00", "#56B4E9", "#009E73", "darkorchid",
                "mistyrose1", "sienna1")




###############################
###############################
######### PhyloSeq ############
###############################
### Prepare count matrix. Filter taxa that have read count < 5 reads 
countTable = read.csv(file="../taxProfiles/CF12/LungTxCF12_countMatrix.csv", header=TRUE, check.names = FALSE, row.names = 1)
tmp = as.matrix(countTable)
# row.names(tmp)=tmp[,1]
# tmp=tmp[,-1]
countTable_2 = apply(tmp, 2, function(x){
  #x[(x/sum(x)) < 0.01] = 0 
  x[x < 5] = 0 
  x
})
countTable_3 = countTable_2[rowSums(countTable_2)!=0,]

## Prepare Taxa matrix
taxaTable = read.csv(file="../taxProfiles/CF12/LungTxCF12_OTUMatrix.csv", header=TRUE, check.names = FALSE, row.names = 1)
taxaTable_2 = as.matrix(taxaTable)

# Prepare metadata matrix
meta = read.csv(file="metadata/CF_12_v2.4.csv", header = TRUE)

## Prepare Phyloseq data object
OTU = otu_table(countTable_3, taxa_are_rows = TRUE)
TAX = tax_table(taxaTable_2)
META = sample_data(meta)
sample_names(META) = meta$SampleIDFP
physeq = phyloseq(OTU, TAX, META)
physeq = subset_taxa(physeq, Superkingdom != "")


## Prepare unfiltered Phyloseq data object
OTU_unfilt = otu_table(tmp, taxa_are_rows = TRUE)
physeq_unfilt = phyloseq(OTU_unfilt, TAX, META)
physeq_unfilt = subset_taxa(physeq_unfilt, Superkingdom != "")


#############################################
### Visualize the 4 superkingdoms proprtion
#############################################
physeqRA = transform_sample_counts(physeq, function(x) x / sum(x))
physeqRA.glom = tax_glom(physeqRA, "Superkingdom")

jpeg("LgTxCF12_all_superkingdom_RA.jpg", res = 300, height = 15, width = 30, units = 'cm')
p = plot_bar(physeqRA.glom, x = "id", fill = "Superkingdom")#, facet_grid = ~Outcome)
p + geom_bar(stat="identity", position="stack") + theme(legend.position="bottom") + 
  ggtitle("Superkingdoms' Proportions") + 
  labs(y = "Relative Abundance", x = "Sample") +
  theme(axis.text.x = element_text(colour="black", size=7, angle=45, hjust=1, vjust=1, face="plain"),
        axis.text.y = element_text(colour="black", size=10, angle=0, hjust=0.5, vjust=0.5, face="plain"),
        axis.title.x = element_text(colour="black", size=12, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=12, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=12, face="bold"), legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
dev.off()



####################################################
## Extract all microbial & Bacteria & Virus species 
####################################################
# filtered
eukaryota = subset_taxa(physeq, Superkingdom == "Eukaryota")
microbial = subset_taxa(physeq, Superkingdom != "Eukaryota")
bacteria = subset_taxa(microbial, Superkingdom == "Bacteria")
virus = subset_taxa(microbial, Superkingdom == "Viruses")
archaea = subset_taxa(microbial, Superkingdom = "Archaea")

## Unfiltered
eukaryota_unfilt = subset_taxa(physeq_unfilt, Superkingdom == "Eukaryota")
microbial_unfilt = subset_taxa(physeq_unfilt, Superkingdom != "Eukaryota")
bacteria_unfilt = subset_taxa(microbial_unfilt, Superkingdom == "Bacteria")
virus_unfilt = subset_taxa(microbial_unfilt, Superkingdom == "Viruses")
archaea_unfilt = subset_taxa(microbial_unfilt, Superkingdom = "Archaea")

## Extract Human reads
human = subset_taxa(physeq, Species == "Homo_sapiens")
human_df = as.data.frame(apply(otu_table(human), 2, sum))
colnames(human_df)[1]="count"
write.csv(human_df, file="Number of Human reads per sample.csv")

## Some numbers analysis on bacteria

View(tax_table(bacteria))

length(unique(tax_table(bacteria)[,"Phylum"]))
length(unique(tax_table(bacteria)[,"Family"]))
length(unique(tax_table(bacteria)[,"Genus"]))
length(unique(tax_table(bacteria)[,"Species"]))


########################################
###### Rarefaction Curve Function ######
########################################
calculate_rarefaction_curves <- function(psdata, measures, depths) {
  require('plyr') # ldply
  require('reshape2') # melt
  
  estimate_rarified_richness <- function(psdata, measures, depth) {
    if(max(sample_sums(psdata)) < depth) return()
    psdata <- prune_samples(sample_sums(psdata) >= depth, psdata)
    
    rarified_psdata <- rarefy_even_depth(psdata, depth, verbose = FALSE)
    alpha_diversity <- estimate_richness(rarified_psdata, measures = measures)
    
    # as.matrix forces the use of melt.array, which includes the Sample names (rownames)
    molten_alpha_diversity <- melt(as.matrix(alpha_diversity), varnames = c('Sample', 'Measure'), value.name = 'Alpha_diversity')
    molten_alpha_diversity
  }
  
  names(depths) <- depths # this enables automatic addition of the Depth to the output by ldply
  rarefaction_curve_data <- ldply(depths, estimate_rarified_richness, psdata = psdata, measures = measures, .id = 'Depth', .progress = ifelse(interactive(), 'text', 'none'))
  
  # convert Depth from factor to numeric
  rarefaction_curve_data$Depth <- as.numeric(levels(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]
  rarefaction_curve_data
}


### Generate Rarefied Data ##list of richness measures ###Depths to rarefy at 
#### how many times to rarefy at each depth
rarefaction_curve_data = calculate_rarefaction_curves(microbial, c('Observed'), rep(c(1, 1000, 2500, 3750, 5000, 7500, 10000, 15000), each = 10))
summary(rarefaction_curve_data)
max(sample_sums(microbial))
rarefaction_curve_data_summary = ddply(rarefaction_curve_data, c('Depth', 'Sample', 'Measure'), 
                                        summarise, Alpha_diversity_mean = mean(Alpha_diversity), 
                                        Alpha_diversity_sd = sd(Alpha_diversity))
#colnames(rarefaction_curve_data_summary)[2] = "SampleIDFP"
rarefaction_curve_data_summary_verbose = merge(rarefaction_curve_data_summary, data.frame(sample_data(microbial)), by.x = 'Sample', by.y = 'row.names')


jpeg("LgTxCF12_RarefactionCurve_0.4.jpg", res = 300, height = 15, width = 20, units = 'cm')
ggplot( data = rarefaction_curve_data_summary_verbose,
             mapping = aes(x = Depth, y = Alpha_diversity_mean, 
                ymin = Alpha_diversity_mean - Alpha_diversity_sd,
                ymax = Alpha_diversity_mean + Alpha_diversity_sd, colour = AnnotatedID,
                group = Sample)) + 
  geom_line() + geom_pointrange(size = 0.5, alpha = 0.5) + 
  scale_color_manual(values = cbbPalette_12) +
  facet_wrap(facets = ~ Measure,  scales = 'free_y') +
  theme_bw() +
  theme(axis.text.x = element_text(colour="black", size=10, angle=0, hjust=0.5, vjust=0, face="plain"),
        axis.text.y = element_text(colour="black", size=10, angle=0, hjust=0.5, vjust=0.5, face="plain"),
        axis.title.x = element_text(colour="black", size=12, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=12, angle=90, hjust=.5, vjust=.5, face="bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position="top",
        strip.text = element_text(size=12, face = "bold"), 
        legend.text = element_text(size=14, face="bold"), legend.title = element_blank()) +
        guides(colour = guide_legend(nrow = 2))
dev.off()



############################################
##### Alpha Diversity using phyloseq #######
############################################
## filtered taxa diveresity
plot_richness(physeq, x = "Outcome", color = "Outcome", measures = c("shannon", "InvSimpson", "Fisher")) + geom_boxplot() + geom_jitter()
plot_richness(microbial, x = "Outcome", color = "Outcome", measures = c("shannon", "InvSimpson", "Fisher")) + geom_boxplot() + geom_jitter()

## unfiltered taxa diversity: microbial unfiltered works better
plot_richness(physeq_unfilt, x = "Outcome", color = "Outcome", measures = c("shannon", "InvSimpson", "Fisher")) + geom_boxplot() + geom_jitter()
plot_richness(microbial_unfilt, x = "Outcome", color = "Outcome", measures = c("shannon", "InvSimpson", "Fisher")) + geom_boxplot() + geom_jitter()


jpeg("LgTxCF12_alpha_diversity_violin_unfiltered_fisher.jpg", res = 300, height = 7, width = 5, units = 'cm')
plot_richness(microbial_unfilt, x = "Outcome", color = "Outcome", 
              measures = c("Fisher")) + 
  geom_violin(aes(fill = Outcome), trim = FALSE) + 
  #geom_boxplot() +
  # geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, position=position_dodge(1)) + 
  ggtitle("")+#"Fisher Diversity Index per Group") +  
  scale_fill_manual(values=c( "red", "midnightblue"),  breaks=c("bos", "non-bos")) + 
  #scale_color_manual(values = cbPalette) +
  #scale_fill_brewer(palette = "Set3")
  #labs(y = "Relative Abundance", x = "Time Point") + 
  theme_bw() + stat_summary(fun.y=mean, geom="point", size=2, color="black") + geom_boxplot(width=0.1) +
  #scale_y_continuous(breaks = c(0, 1)) + 
  theme(axis.text.x = element_text(colour="black", size=10, angle=0, hjust=0.5, vjust=0, face="bold"),
        axis.text.y = element_text(colour="black", size=10, angle=0, hjust=0.5, vjust=0.5, face="bold"),
        axis.title.x = element_text(colour="black", size=10, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=10, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=8, face="plain"), 
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position="none",
        strip.text = element_text(size=12, face = "bold"))
dev.off()


######################################
######### Visualize diversity per subject
######################################
#sub2 = subset_samples(microbial_unfilt, Outcome == "bos")
jpeg("LgTxCF12_alpha_diversity_perSubject_Fisher.jpg", res = 300, height = 10, width = 10, units = 'cm')
plot_richness(microbial_unfilt, x = "AnnotatedID", color = "Outcome", 
              measures = c("Fisher")) + 
  #geom_violin(aes(fill = TimePoint), trim = FALSE) + 
  #geom_boxplot() +
  # geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, position=position_dodge(1)) + 
  ggtitle("") +  
  scale_fill_manual(values=c( "red", "midnightblue"),  breaks=c("bos", "non-bos")) + 
  scale_color_manual(values=c( "red", "midnightblue"),  breaks=c("bos", "non-bos")) +
  #scale_fill_brewer(palette = "Set3")
  labs( x = "Subject") + 
  theme_bw() + stat_summary(fun.y=mean, geom="point", size=1, color="black") + geom_boxplot(width=0.1) +
  #scale_y_continuous(breaks = c(0, 1)) + 
  theme(axis.text.x = element_text(colour="black", size=10, angle=45, hjust=1, vjust=1, face="bold"),
        axis.text.y = element_text(colour="black", size=10, angle=0, hjust=0.5, vjust=0.5, face="bold"),
        axis.title.x = element_text(colour="black", size=10, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=10, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=15, face="plain"), 
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position="none",
        strip.text = element_text(size=12, face = "bold")) +
  guides(colour = guide_legend(nrow = 1))
dev.off()


############################################################################
######### Visualize diversity Tragectories per subject  ################
############################################################################
### Plot for the first two timepoints and last timepoint
## merge meta and richness
richness = estimate_richness(microbial_unfilt)
richness$SampleIDFP = rownames(richness)
meta_diversity = merge(meta,richness, by="SampleIDFP")
View(meta_diversity)
#write.csv(meta_diversity, file="meta_diversity.csv", row.names = FALSE)


######################
######################
# Get median of of diversity

meta_diversity$Fisher
meta_diversity$Fisher_Diversity

myFun <- function(x) {
  median = median(x)
}

x = tapply(meta_diversity$Fisher_Diversity, meta_diversity$SUBJECT.WUTX, myFun)
x = as.data.frame(x)

x$SUBJECT.WUTX = rownames(x)
y = unique(meta_diversity[,c("SUBJECT.WUTX", "OutcomeMod")])
z = merge(x,y, by = "SUBJECT.WUTX")
View(x)
View(y)
View(z)

wilcox.test(x~OutcomeMod, z)



median(z[z[,"OutcomeMod"]=="BO",]$x)
median(z[z[,"OutcomeMod"]=="non-BO",]$x)



## BO
sub2 = subset(meta_diversity, Outcome == "BO")
sub3 = subset(sub2, bin_status_mod == "A50" | bin_status_mod == "B100" | bin_status == "last")
jpeg("LgTxCF12_alpha_diversity_unfiltered_BOS_3timepoints.jpg", res = 300, height = 10, width = 10, units = 'cm')
ggplot(sub3, aes(x = bin_status_mod, y=InvSimpson, color = SUBJECT.WUTX, group = SUBJECT.WUTX)) + 
  geom_point(size=4, alpha=0.8) + geom_line(alpha=0.5, size=2) +
  ggtitle("Diversity Index for BOS") +  
  scale_fill_manual(values = cbbPalette_12) +
  scale_color_manual(values=cbbPalette_12) + 
  labs(x = "Timepoint") + 
  theme_bw() +  
  theme(axis.text.x = element_text(colour="black", size=10, angle=0, hjust=0.5, vjust=0, face="bold"),
        axis.text.y = element_text(colour="black", size=10, angle=0, hjust=0.5, vjust=0.5, face="bold"),
        axis.title.x = element_text(colour="black", size=12, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=12, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=10, face="plain"), 
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
        legend.position="right",
        strip.text = element_text(size=12, face = "bold")) +
  guides(colour = guide_legend(ncol = 1))
dev.off()




## Nonbos
sub2 = subset(meta_diversity, Outcome == "non-bos")
sub3 = subset(sub2, bin_status_mod == "A50" | bin_status_mod == "B100" | bin_status == "last")
jpeg("LgTxCF12_alpha_diversity_unfiltered_nonBOS_3timepoints.jpg", res = 300, height = 10, width = 10, units = 'cm')
ggplot(sub3, aes(x = bin_status_mod, y=InvSimpson, color = SUBJECT.WUTX, group = SUBJECT.WUTX)) + 
  geom_point(size=4, alpha=0.8) + geom_line(alpha=0.5, size=2) +  
  ggtitle("Diversity Index for non-BOS") +  
  scale_fill_manual(values = cbbPalette_12) +
  scale_color_manual(values=cbbPalette_12) + 
  labs(x = "Timepoint") + 
  theme_bw() +  
  theme(axis.text.x = element_text(colour="black", size=10, angle=0, hjust=0.5, vjust=0, face="bold"),
        axis.text.y = element_text(colour="black", size=10, angle=0, hjust=0.5, vjust=0.5, face="bold"),
        axis.title.x = element_text(colour="black", size=12, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=12, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=10, face="plain"), 
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
        legend.position="right",
        strip.text = element_text(size=12, face = "bold")) +
  guides(colour = guide_legend(ncol = 1))
dev.off()










################ Using Facet
#sub2 = subset(meta_diversity, Outcome == "BO")
sub3 = subset(meta, bin_status_mod == "A50" | bin_status_mod == "B100" | bin_status == "last")
jpeg("LgTxCF12_alpha_diversity_unfiltered_BOS_3timepoints_bos_nonbos.jpg", res = 300, height = 10, width = 20, units = 'cm')
ggplot(sub3, aes(x = bin_status_mod, y=Fisher_Diversity, color = AnnotatedID, group = AnnotatedID)) + 
  geom_point(size=4, alpha=0.8) + geom_line(alpha=0.5, size=2) +
  ggtitle("") +  
  scale_fill_manual(values = cbbPalette_12) +
  scale_color_manual(values=cbbPalette_12) + 
  labs(x = "Timepoint") + 
  theme_bw() +  
  theme(axis.text.x = element_text(colour="black", size=10, angle=0, hjust=0.5, vjust=0, face="bold"),
        axis.text.y = element_text(colour="black", size=10, angle=0, hjust=0.5, vjust=0.5, face="bold"),
        axis.title.x = element_text(colour="black", size=12, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=12, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=10, face="plain"), 
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
        legend.position="right",
        strip.text = element_text(size=12, face = "bold")) +
  guides(colour = guide_legend(ncol = 1)) + facet_grid(. ~ OutcomeMod)
dev.off()










### For all of them 
# bos
sub2 = subset_samples(microbial_unfilt, Outcome == "BO")
jpeg("LgTxCF12_alpha_diversity_unfiltered_BOS_trajectory.jpg", res = 300, height = 15, width = 20, units = 'cm')
plot_richness(sub2, x = "TimeAfterTransplantation", color = "AnnotatedID", 
              measures = c("Fisher")) + geom_point(size=6, alpha=0.5) + geom_line(alpha=0.5, size=2, group = "AnnotatedID") +
  #geom_violin(aes(fill = TimePoint_Status), trim = FALSE) + 
  #geom_boxplot() +
  # geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, position=position_dodge(1)) + 
  ggtitle("Diversity") +  
  scale_fill_manual(values = cbbPalette_12) +
  scale_color_manual(values=cbbPalette_12) + 
  #scale_color_manual(values = cbbPalette_12) +
  #scale_fill_brewer(palette = "Set3")
  labs(x = "Timepoint") + 
  theme_bw() + #stat_summary(fun.y=mean, geom="point", size=1, color="black") + #geom_boxplot(width=0.1) +
  #scale_y_continuous(breaks = c(0, 1)) + 
  theme(axis.text.x = element_text(colour="black", size=10, angle=0, hjust=0.5, vjust=0, face="bold"),
        axis.text.y = element_text(colour="black", size=10, angle=0, hjust=0.5, vjust=0.5, face="bold"),
        axis.title.x = element_text(colour="black", size=12, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=12, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=10, face="plain"), 
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
        legend.position="right",
        strip.text = element_text(size=12, face = "bold")) +
  guides(colour = guide_legend(ncol = 1)) + facet_grid(AnnotatedID ~ Outcome)#facet_grid(SUBJECT.WUTX ~ .)
dev.off()


## nonbos
sub2 = subset_samples(microbial_unfilt, Outcome == "non-BO")
jpeg("LgTxCF12_alpha_diversity_unfiltered_nonBOS_trajectory.jpg", res = 300, height = 20, width = 20, units = 'cm')
plot_richness(sub2, x = "TimeAfterTransplantation", color = "AnnotatedID", 
              measures = c("Fisher")) + geom_point(size=6, alpha=0.5) + geom_line(alpha=0.5, size=2, group = "AnnotatedID") +
  #geom_violin(aes(fill = TimePoint_Status), trim = FALSE) + 
  #geom_boxplot() +
  # geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, position=position_dodge(1)) + 
  ggtitle("Diversity") +  
  scale_fill_manual(values = cbbPalette_12) +
  scale_color_manual(values=cbbPalette_12) + 
  #scale_color_manual(values = cbbPalette_12) +
  #scale_fill_brewer(palette = "Set3")
  labs(x = "Timepoint") + 
  theme_bw() + #stat_summary(fun.y=mean, geom="point", size=1, color="black") + #geom_boxplot(width=0.1) +
  #scale_y_continuous(breaks = c(0, 1)) + 
  theme(axis.text.x = element_text(colour="black", size=10, angle=0, hjust=0.5, vjust=0, face="bold"),
        axis.text.y = element_text(colour="black", size=10, angle=0, hjust=0.5, vjust=0.5, face="bold"),
        axis.title.x = element_text(colour="black", size=12, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=12, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=10, face="plain"), 
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
        legend.position="right",
        strip.text = element_text(size=12, face = "bold")) +
  guides(colour = guide_legend(ncol = 1)) + facet_grid(AnnotatedID ~ Outcome)#facet_grid(SUBJECT.WUTX ~ .)
dev.off()


#######################################
### Statistics on diversity
#######################################
### TODO: Add pvalue to teh graphs
richness = estimate_richness(microbial_unfilt)
mannwhitney_test = t(sapply(richness, function(x) unlist(wilcox.test(x ~ sample_data(microbial_unfilt)$Outcome) [c("estimate","p.value","statistic","conf.int")])))
t(sapply(richness, function(x) unlist(wilcox.test(x ~ sample_data(microbial_unfilt)$Outcome))))

# Shannon    "493"       "0.0042651927402922" "0"                       "two.sided"
# Simpson    "533"       "0.013186185727668"  "0"                       "two.sided"
# InvSimpson "533"       "0.013186185727668"  "0"                       "two.sided"
# Fisher     "628"       "0.114515817365666"  "0"                       "two.sided"







#############################
######## Do MetaLonDA #######
#############################
install.packages("MetaLonDA")
library(devtools)
install_github("aametwally/MetaLonDA", ref = "v1.1.2")
library(MetaLonDA)
library(zoo)


## MetaLonDA on Diversity
count = meta_diversity$Fisher
time = meta_diversity$TimeAfterTransplantation
id = meta_diversity$AnnotatedID
group = meta_diversity$Outcome
points = seq(1, 1025, length.out = 100)
output.metalonda.f1 = metalonda(Count = as.matrix(count), Time = time, Group = group,
                                ID = id, n.perm = 1000, fit.method = "lowess", points = points,
                                text = "diversity", parall = FALSE, pvalue.threshold = 0.1,     
                                adjust.method = "BH", col = c("firebrick", "blue"), ylabel = "Diversity")




## MetaLonDA on FEV
count = meta_diversity$FEV1_Per_Predicted
time = meta_diversity$TimeAfterTransplantation
id = meta_diversity$AnnotatedID
group = meta_diversity$Outcome
points = seq(1, 1025, length.out = 100)
output.metalonda.f1 = metalonda(Count = as.matrix(count), Time = time, Group = group,
                                ID = id, n.perm = 1000, fit.method = "lowess", points = points,
                                text = "FEV1 Percent Predicted", parall = FALSE, pvalue.threshold = 0.1,     
                                adjust.method = "BH", col = c("firebrick", "blue"), ylabel = "FEV1 (%) Predicted")



## MetaLonDA on FVC
count = meta_diversity$FVC_Per_Predicted
time = meta_diversity$TimeAfterTransplantation
id = meta_diversity$AnnotatedID
group = meta_diversity$Outcome
points = seq(1, 1025, length.out = 100)
output.metalonda.f1 = metalonda(Count = as.matrix(count), Time = time, Group = group,
                                ID = id, n.perm = 1000, fit.method = "lowess", points = points,
                                text = "FVC Percent Predicted", parall = FALSE, pvalue.threshold = 0.1,     
                                adjust.method = "BH", col = c("firebrick", "blue"), ylabel = "FVC (%) Predicted")


## MetaLonDA on FEV1/FVC
count = meta_diversity$FEV1_FVC_Ratio
time = meta_diversity$TimeAfterTransplantation
id = meta_diversity$AnnotatedID
group = meta_diversity$Outcome
points = seq(1, 1025, length.out = 100)
output.metalonda.f1 = metalonda(Count = as.matrix(count), Time = time, Group = group,
                                ID = id, n.perm = 1000, fit.method = "lowess", points = points,
                                text = "FEV1 over FVC", parall = FALSE, pvalue.threshold = 0.1,     
                                adjust.method = "BH", col = c("firebrick", "blue"), ylabel = "FEV1/FVC")




## MetaLonDA on LavageTotalCellCountcells_per_mcL
count = meta_diversity$LavageTotalCellCountcells_per_mcL
time = meta_diversity$TimeAfterTransplantation
id = meta_diversity$AnnotatedID
group = meta_diversity$Outcome
points = seq(1, 460, length.out = 100)
output.metalonda.f1 = metalonda(Count = as.matrix(count), Time = time, Group = group,
                                ID = id, n.perm = 1000, fit.method = "lowess", points = points,
                                text = "Lavage Total CellCount cells_per_mcL", parall = FALSE, pvalue.threshold = 0.1,     
                                adjust.method = "BH", col = c("firebrick", "blue"), ylabel = "cells_per_mcL")


## MetaLonDA on Neutrophil_perc
count = meta_diversity$Neutrophil_perc
time = meta_diversity$TimeAfterTransplantation
id = meta_diversity$AnnotatedID
group = meta_diversity$Outcome
points = seq(1, 1025, length.out = 100)
output.metalonda.f1 = metalonda(Count = as.matrix(count), Time = time, Group = group,
                                ID = id, n.perm = 1000, fit.method = "lowess", points = points,
                                text = "Neutrophil_perc", parall = FALSE, pvalue.threshold = 0.1,     
                                adjust.method = "BH", col = c("firebrick", "blue"), ylabel = "Neutrophil_perc")



## MetaLonDA on Lymphocyte_perc
count = meta_diversity$Lymphocyte_perc
time = meta_diversity$TimeAfterTransplantation
id = meta_diversity$AnnotatedID
group = meta_diversity$Outcome
points = seq(1, 1025, length.out = 100)
output.metalonda.f1 = metalonda(Count = as.matrix(count), Time = time, Group = group,
                                ID = id, n.perm = 1000, fit.method = "lowess", points = points,
                                text = "Lymphocyte_perc", parall = FALSE, pvalue.threshold = 0.1,     
                                adjust.method = "BH", col = c("firebrick", "blue"), ylabel = "Lymphocyte_perc")





## MetaLonDA on RecentTacrolimusTroughLevel_ng_per_mL
count = meta_diversity$RecentTacrolimusTroughLevel_ng_per_mL
time = meta_diversity$TimeAfterTransplantation
id = meta_diversity$AnnotatedID
group = meta_diversity$Outcome
points = seq(1, 1025, length.out = 100)
output.metalonda.f1 = metalonda(Count = as.matrix(count), Time = time, Group = group,
                                ID = id, n.perm = 1000, fit.method = "lowess", points = points,
                                text = "RecentTacrolimusTroughLevel_ng_per_mL", parall = FALSE, pvalue.threshold = 0.1,     
                                adjust.method = "BH", col = c("firebrick", "blue"), ylabel = "RecentTacrolimusTroughLevel_ng_per_mL")



## MetaLonDA on RecentMycophenolicAcidLevel_mcg_per_mL
count = meta_diversity$RecentMycophenolicAcidLevel_mcg_per_mL
time = meta_diversity$TimeAfterTransplantation
id = meta_diversity$AnnotatedID
group = meta_diversity$Outcome
points = seq(1, 1025, length.out = 100)
output.metalonda.f1 = metalonda(Count = as.matrix(count), Time = time, Group = group,
                                ID = id, n.perm = 1000, fit.method = "lowess", points = points,
                                text = "RecentMycophenolicAcidLevel_mcg_per_mL", parall = FALSE, pvalue.threshold = 0.1,     
                                adjust.method = "BH", col = c("firebrick", "blue"), ylabel = "RecentMycophenolicAcidLevel_mcg_per_mL")






####################################################
#### Metalonda on new file from Cristian/Ferkol ####
####################################################
bal_acr = read.csv("BALcells_ACR_Metadata_CF_v1.5.csv", header = TRUE)

## LavageTotalCellCountcells_per_mcL
count = bal_acr$LavageTotalCellCountcells_per_mcL
time = bal_acr$TimeAfterTransplantation
id = bal_acr$AnnotatedID
group = bal_acr$Outcome
points = seq(1, 460, length.out = 100)
output.metalonda.f1 = metalonda(Count = as.matrix(count), Time = time, Group = group,
                                ID = id, n.perm = 1000, fit.method = "lowess", points = points,
                                text = "LavageTotalCellCountcells_per_mcL", parall = FALSE, pvalue.threshold = 0.05,     
                                adjust.method = "BH", col = c("firebrick", "blue"), ylabel = "cells_per_mcL")


## Neutrophil_perc
count = bal_acr$Neutrophil_perc
time = bal_acr$TimeAfterTransplantation
id = bal_acr$AnnotatedID
group = bal_acr$Outcome
points = seq(1, 460, length.out = 100)
output.metalonda.f1 = metalonda(Count = as.matrix(count), Time = time, Group = group,
                                ID = id, n.perm = 1000, fit.method = "lowess", points = points,
                                text = "Neutrophil_perc", parall = FALSE, pvalue.threshold = 0.05,     
                                adjust.method = "BH", col = c("firebrick", "blue"), ylabel = "Neutrophil_perc")



## Lymphocyte_perc
count = bal_acr$Lymphocyte_perc
time = bal_acr$TimeAfterTransplantation
id = bal_acr$AnnotatedID
group = bal_acr$Outcome
points = seq(1, 460, length.out = 100)
output.metalonda.f1 = metalonda(Count = as.matrix(count), Time = time, Group = group,
                                ID = id, n.perm = 1000, fit.method = "lowess", points = points,
                                text = "Lymphocyte_perc", parall = FALSE, pvalue.threshold = 0.05,     
                                adjust.method = "BH", col = c("firebrick", "blue"), ylabel = "Lymphocyte_perc")







###########################################
##### Correlate/Linear regression Diversity with Neutrophil_perc
###########################################
fit_bos <- lm(Neutrophil_perc ~ Fisher, subset = Outcome == "BO",data = meta_diversity)
summary(fit_bos)


fit_non_bos <- lm(Neutrophil_perc ~ Fisher, subset = Outcome == "non-BO",data = meta_diversity)
summary(fit_non_bos)




# p-value = 3.54e-05
# Multiple R-squared:  0.2027,	Adjusted R-squared:  0.1922


jpeg("LgTxCF12_Neutrophils__diversity_pergroup.jpg", res = 300, height = 8, width = 15, units = 'cm')
ggplot(meta_diversity, aes(x = Fisher, y = Neutrophil_perc, color = Outcome)) + 
  geom_point() +
  stat_smooth(method = "lm", aes(color=Outcome)) +
  ggtitle("Neutrophils % vs Diversity") +  
  scale_fill_manual(values=c( "red", "midnightblue"),  breaks=c("BO", "non-BO")) + 
  scale_color_manual(values=c( "red", "midnightblue"),  breaks=c("BO", "non-BO")) +
  #scale_fill_manual(values = cbbPalette_12) +
  #scale_color_manual(values=cbbPalette_12) + 
  #scale_color_manual(values = cbbPalette_12) +
  #scale_fill_brewer(palette = "Set3")
  labs(x = "Diversity", y = "Neutrophils %") + 
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






####### Lymphocytes and Diversity
fit_bos <- lm(Lymphocyte_perc ~ Fisher, subset = Outcome == "BO",data = meta_diversity)
summary(fit_bos)


fit_non_bos <- lm(Lymphocyte_perc ~ Fisher, subset = Outcome == "non-BO",data = meta_diversity)
summary(fit_non_bos)


jpeg("LgTxCF12_Lymphocytes__diversity_pergroup.jpg", res = 300, height = 8, width = 15, units = 'cm')
ggplot(meta_diversity, aes(x = Fisher, y = Lymphocyte_perc, color = Outcome)) + 
  geom_point() +
  stat_smooth(method = "lm", aes(color=Outcome)) +
  ggtitle("Lymphocytes % vs Diversity") +  
  scale_fill_manual(values=c( "red", "midnightblue"),  breaks=c("BO", "non-BO")) + 
  scale_color_manual(values=c( "red", "midnightblue"),  breaks=c("BO", "non-BO")) +
  #scale_fill_manual(values = cbbPalette_12) +
  #scale_color_manual(values=cbbPalette_12) + 
  #scale_color_manual(values = cbbPalette_12) +
  #scale_fill_brewer(palette = "Set3")
  labs(x = "Diversity", y = "Lymphocytes %") + 
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



####### Nutrophils and Diversity (categorized by microbiology status)
meta_diversity$culture ="+"
meta_diversity$culture[which(meta_diversity$LavageMicrobiology == "NoSignificantGrowth" | meta_diversity$LavageMicrobiology == "No significant growth")] = "-"


fit_negCult <- lm(Neutrophil_perc ~ Fisher, subset = culture == "-",data = meta_diversity)
summary(fit_negCult)


fit_posCult <- lm(Neutrophil_perc ~ Fisher, subset = culture == "+",data = meta_diversity)
summary(fit_posCult)




jpeg("LgTxCF12_Neutrophils__diversity_perCulture.jpg", res = 300, height = 8, width = 15, units = 'cm')
ggplot(meta_diversity, aes(x = Fisher, y = Neutrophil_perc, color = culture)) + 
  geom_point() +
  stat_smooth(method = "lm", aes(color=culture)) +
  ggtitle("Neutrophils % vs Diversity") +  
  scale_fill_manual(values=c( "darkgreen", "orange"),  breaks=c("+", "-")) + 
  scale_color_manual(values=c("darkgreen", "orange"),  breaks=c("+", "-")) +
  #scale_fill_manual(values = cbbPalette_12) +
  #scale_color_manual(values=cbbPalette_12) + 
  #scale_color_manual(values = cbbPalette_12) +
  #scale_fill_brewer(palette = "Set3")
  labs(x = "Diversity", y = "Neutrophils %") + 
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
  guides(colour = guide_legend(ncol = 1)) + facet_grid(. ~ culture)
dev.off()





# plot(Neutrophil_perc ~ Fisher, data = meta_diversity)
# abline(fit1)



# ggplot(meta_diversity, aes(Fisher, Neutrophil_perc, color=Outcome)) +
#   stat_smooth(method="loess", formula=y~x,
#               alpha=0.2, size=2, aes(fill=Outcome)) +
#   geom_point(position=position_jitter(height=0.03, width=0)) +
#   xlab("Diversity") + ylab("Neutrophil %")




##############################################
##### Logistic regression Outcome ~ diversity
##############################################
model <- glm(Outcome ~ InvSimpson+Fisher,family=binomial(link='logit'),data=meta_diversity)
summary(model)
anova(model, test="Chisq")

## Assess predictive power of the model 
fitted.results <- predict(model,newdata=meta_diversity,type='response')
fitted.results <- ifelse(fitted.results > 0.5, "bos", "non-bos")

misClasificError <- mean(fitted.results != meta_diversity$Outcome)
print(paste('Accuracy',1-misClasificError))

## plot ROC
install.packages("ROCR")
library(ROCR)

fitted.results <- predict(model,newdata=meta_diversity,type='response')
pr <- prediction(fitted.results, meta_diversity$Outcome)
prf <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf)


auc <- performance(pr, measure = "auc")
auc <- auc@y.values[[1]]
print(paste('AUC = ', auc))






#### COnvert BOS and non-BOS to nuemric 1, 0:
meta_diversity$binaryOutcome = 1
meta_diversity$binaryOutcome[which(meta_diversity$Outcome=="non-bos")] = 0
#write.csv(meta_diversity, file="meta_diversity.csv", row.names = FALSE)


## plot Lowess
ggplot(meta_diversity, aes(Fisher, binaryOutcome, color=RecipientGender)) +
  stat_smooth(method="loess", formula=y~x,
              alpha=0.2, size=2, aes(fill=RecipientGender)) +
  geom_point(position=position_jitter(height=0.03, width=0)) +
  xlab("Diversity") + ylab("Pr(BOS)")


## plot GLM
ggplot(meta_diversity, aes(Fisher, binaryOutcome, color=RecipientGender)) +
  stat_smooth(method="glm", family=binomial, formula=y~x,
              alpha=0.2, size=2, aes(fill=RecipientGender)) +
  geom_point(position=position_jitter(height=0.03, width=0)) +
  xlab("Diversity") + ylab("Pr(BOS)")



## GLM Model
model <- glm(binaryOutcome ~ (InvSimpson + Fisher + Shannon)^3,family=binomial(link='logit'),data=meta_diversity)
summary(model)
anova(model, test="Chisq")
contrasts(meta_diversity$Outcome)


## http://ddar.datavis.ca/pages/extra/titanic-glm-ex.pdf
#install.packages("effects")

# titanic.eff2 <- allEffects(model)
# names(titanic.eff2)
# 
# titanic.eff2a <- allEffects(model,
#                             typical=median)#,
#                             #given.values=c(InvSimpson=1/3, Fisher=1/3, Shannon=0.5)
# #)
# 
# #plot(model, ask=TRUE)
# ticks <- list(at=c(.01, .05, seq(.1, .9, by=.2), .95, .99))
# plot(titanic.eff2[1], ticks=ticks, multiline=TRUE, ci.style="bars"), key=list(x=.7, y=.95))
# 


###########################################
## Survival analysis related to Diversity 
###########################################
library(survminer)
library(survival)

meta_diversity = read.csv("meta_diversity.csv")
subset_survival = unique(meta_diversity[, c("SUBJECT.WUTX", "Outcome", "SurvivalTime", "RecipientGender")])
subset_survival$SurvObj =  with(subset_survival, Surv(SurvivalTime, Outcome == "bos"))

fitsex <- survfit(SurvObj ~ RecipientGender, data = subset_survival, conf.type = "log-log")
plot(fitsex)


jpeg("LgTxCF12_KaplanMeir_Gender.jpg", res = 300, height = 20, width = 20, units = 'cm')
ggsurvplot(fitsex, data = subset_survival, risk.table = TRUE, conf.int = TRUE, pval = TRUE)
dev.off()


fit_all <- survfit(SurvObj ~ 1, data = subset_survival)
ggsurvplot(fit_all, data = subset_survival, risk.table = TRUE, conf.type = "log-log")




#######################################
### Subset the top 5 bacteria phyla 
#######################################
### TODO: Needs to agglomorate
# top5ph = sort(tapply(taxa_sums(bacteria), tax_table(bacteria)[, "Phylum"], sum), decreasing = TRUE)[1:5]
# #names(top5ph)[2] = "virus"
# microbial.top5.phylum = subset_taxa(microbial, Phylum %in% names(top5ph))
# microbial.top5.phylum.RA = transform_sample_counts(microbial.top5.phylum, function(x) x / sum(x))



## optimize
top5ph = sort(tapply(taxa_sums(bacteria), tax_table(bacteria)[, "Phylum"], sum), decreasing = TRUE)[1:5]
top4ph = top5ph[-4]
bacteria.RA = transform_sample_counts(bacteria, function(x) x / sum(x))
bacteria.top4.phylum.RA = subset_taxa(bacteria.RA, Phylum %in% names(top4ph))


jpeg("LUTX_top_4_bacterial_phyla.jpg", res = 300, height = 20, width = 12, units = 'cm')
p = plot_bar(bacteria.top4.phylum.RA, x = "TimePoint", fill="Phylum", 
             facet_grid= AnnotatedID~.)
p + geom_bar(stat="identity", position="stack") + 
  ggtitle("Top 4 Bacterial Phyla") +  
  scale_fill_manual(values = cbbPalette_9) +
  #scale_fill_brewer(palette = "Set3")
  labs(y = "Relative Abundance", x = "Timepoint") + theme_bw() +
  scale_y_continuous(breaks = c(0, 1)) + 
  theme(axis.text.x = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0, face="bold"),
        axis.text.y = element_text(colour="black", size=8, angle=0, hjust=0.5, vjust=0.5, face="bold"),
        axis.title.x = element_text(colour="black", size=15, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=15, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=8, face="plain"), 
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position="bottom",
        strip.text = element_text(size=6, face = "bold"), strip.background =element_rect(aes(fill=bacteria.top4.phylum.RA$Outcome)))
dev.off()






### Think of doing it as trajectory graph
# jpeg("LUTX_top_4_bacterial_phyla.jpg", res = 300, height = 20, width = 20, units = 'cm')
# p = plot_bar(bacteria.top4.phylum.RA, x = "TimePoint", fill="Phylum", 
#              facet_grid= SUBJECT.WUTX~.)
# p + #geom_bar(stat="identity", position="stack") + 
#   
#   geom_point(size=2, alpha=0.5) + geom_line(alpha=0.5, size=2, group = "SUBJECT.WUTX") +
#   
#   
#   ggtitle("Top 4 Bacterial Phyla") +  
#   scale_fill_manual(values = cbbPalette_9) +
#   #scale_fill_brewer(palette = "Set3")
#   labs(y = "Relative Abundance", x = "Timepoint") + theme_bw() +
#   scale_y_continuous(breaks = c(0, 1)) + 
#   theme(axis.text.x = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0, face="bold"),
#         axis.text.y = element_text(colour="black", size=8, angle=0, hjust=0.5, vjust=0.5, face="bold"),
#         axis.title.x = element_text(colour="black", size=15, angle=0, hjust=.5, vjust=0.5, face="bold"),
#         axis.title.y = element_text(colour="black", size=15, angle=90, hjust=.5, vjust=.5, face="bold"),
#         legend.text=element_text(size=15, face="plain"), 
#         legend.title = element_blank(),
#         plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
#         legend.position="bottom",
#         strip.text = element_text(size=12, face = "bold"))
# dev.off()



################################################
### Subset the top 10 Proteobacteria genera  ###
################################################
proteo = subset_taxa(bacteria, Phylum == "Proteobacteria")
proteo_genus = tax_glom(proteo, taxrank="Genus")
top10genus = sort(tapply(taxa_sums(proteo_genus), tax_table(proteo_genus)[, "Genus"], sum), decreasing = TRUE)[1:10]
top8genus = top10genus[-c(2,8)]
proteo.genus.RA = transform_sample_counts(proteo_genus, function(x) x / sum(x))
proteo.genus.RA.top8 = subset_taxa(proteo.genus.RA, Genus %in% names(top8genus))




jpeg("LgTxCF12_top_8_proteo_genera.jpg", res = 300, height = 20, width = 12, units = 'cm')
p = plot_bar(proteo.genus.RA.top8, x = "TimePoint", fill="Genus", 
             facet_grid= AnnotatedID ~ .)
p + geom_bar(stat="identity", position="stack") + 
  ggtitle("Top 8 Proteobacteria Genera") +  
  scale_fill_manual(values = cbbPalette_10) +
  #scale_fill_brewer(palette = "Set3")
  labs(y = "Relative Abundance", x = "Timepoint") + theme_bw() +
  scale_y_continuous(breaks = c(0, 1)) + 
  theme(axis.text.x = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0, face="bold"),
        axis.text.y = element_text(colour="black", size=8, angle=0, hjust=0.5, vjust=0.5, face="bold"),
        axis.title.x = element_text(colour="black", size=15, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=15, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=8, face="plain"), 
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position="bottom",
        strip.text = element_text(size=6, face = "bold"))
dev.off()





################################################
### Subset the top 10 Firmicutes genera  ###
################################################



firm = subset_taxa(bacteria, Phylum == "Firmicutes")
firm_genus = tax_glom(firm, taxrank="Genus")
top10genus = sort(tapply(taxa_sums(firm_genus), tax_table(firm_genus)[, "Genus"], sum), decreasing = TRUE)[1:10]
top8genus = top10genus[-c(9,10)]
firm.genus.RA = transform_sample_counts(firm_genus, function(x) x / sum(x))
firm.genus.RA.top8 = subset_taxa(firm.genus.RA, Genus %in% names(top8genus))




jpeg("LgTxCF12_top_8_firm_genera.jpg", res = 300, height = 20, width = 12, units = 'cm')
p = plot_bar(firm.genus.RA.top8, x = "TimePoint", fill="Genus", 
             facet_grid= AnnotatedID ~ .)
p + geom_bar(stat="identity", position="stack") + 
  ggtitle("Top 8 Firmicutes Genera") +  
  scale_fill_manual(values = cbbPalette_10) +
  #scale_fill_brewer(palette = "Set3")
  labs(y = "Relative Abundance", x = "Timepoint") + theme_bw() +
  scale_y_continuous(breaks = c(0, 1)) + 
  theme(axis.text.x = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0, face="bold"),
        axis.text.y = element_text(colour="black", size=8, angle=0, hjust=0.5, vjust=0.5, face="bold"),
        axis.title.x = element_text(colour="black", size=15, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=15, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=8, face="plain"), 
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position="bottom",
        strip.text = element_text(size=6, face = "bold"))
dev.off()






################################################
### Subset the top 10 Burkholderia genera  ###
################################################
proteo = subset_taxa(bacteria, Genus == "Burkholderia")
proteo_genus = proteo#tax_glom(proteo, taxrank="Species")
top10genus = sort(tapply(taxa_sums(proteo_genus), tax_table(proteo_genus)[, "Species"], sum), decreasing = TRUE)#[1:10]
top8genus = top10genus# top10genus[-c(2,8)]
proteo.genus.RA = transform_sample_counts(proteo_genus, function(x) x / sum(x))
proteo.genus.RA.top8 = subset_taxa(proteo.genus.RA, Species %in% names(top8genus))



### Save Burk Species

write.csv(otu_table(proteo.genus.RA.top8), file = "Burk_10_Species.csv")

write.csv(tax_table(proteo.genus.RA.top8), file="Burk_10_species_taxonomy.csv")


jpeg("LgTxCF12_Burk_species.jpg", res = 300, height = 20, width = 20, units = 'cm')
p = plot_bar(proteo.genus.RA.top8, x = "TimePoint", fill="Species", 
             facet_grid= AnnotatedID ~ .)
p + geom_bar(stat="identity", position="stack") + 
  ggtitle("Burkholderia Species") +  
  #scale_fill_manual(values = cbbPalette_10) +
  #scale_fill_brewer(palette = "Set3")
  labs(y = "Relative Abundance", x = "Timepoint") + theme_bw() +
  scale_y_continuous(breaks = c(0, 1)) + 
  theme(axis.text.x = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0, face="bold"),
        axis.text.y = element_text(colour="black", size=8, angle=0, hjust=0.5, vjust=0.5, face="bold"),
        axis.title.x = element_text(colour="black", size=15, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=15, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=5, face="plain"), 
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position="right",
        strip.text = element_text(size=6, face = "bold"))
dev.off()




###################################
########  save aggregated 4 phyla for mixed effects model ##############
###################################

top5ph = sort(tapply(taxa_sums(bacteria), tax_table(bacteria)[, "Phylum"], sum), decreasing = TRUE)[1:5]
top4ph = top5ph[-4]
bacteria.top4.phylum.RA = subset_taxa(bacteria, Phylum %in% names(top4ph))
proteo_genus = tax_glom(bacteria.top4.phylum.RA, taxrank="Phylum")
proteo.genus.RA = transform_sample_counts(proteo_genus, function(x) x / sum(x))

write.csv(otu_table(proteo.genus.RA), file = "LgTx_PhylaAbundance.csv")
View(tax_table(proteo.genus.RA))



################## Save Family

families = tax_glom(bacteria, taxrank="Family")
families.RA = transform_sample_counts(families, function(x) x / sum(x))

names = c("31958", "193")
families.RA.subset = subset_taxa(families.RA, Family %in% names(names))


write.csv(otu_table(families.RA), file = "LgTx_FamilyAbundance.csv")
View(tax_table(families.RA))

Pseudonocardiaceae_31958
Rhodospirillaceae_193


f = read.csv("LgTx_FamilyAbundance.csv")
View(f)
meta_family = merge(meta, f, by.x = 'SampleIDFP')
View(meta_family)

write.csv(meta_family, file = "CF_12_v2.4.csv")



# TODO: top 30 genera
top30genera = sort(tapply(taxa_sums(bacteria), tax_table(bacteria)[, "Genus"], sum), decreasing = TRUE)[1:30]



#############################################################
######## plot number of microbial reads per sample ##########
#############################################################
df = as.data.frame(apply(otu_table(microbial), 2, sum))
colnames(df)[1]="count"
write.csv(df, file="Number of Microbial reads per sample.csv")
df$SampleIDFP = rownames(df)



df_merge = merge(df, data.frame(sample_data(microbial)), by.x = 'SampleIDFP', by.y = 'row.names')
df_merge$id_order = paste(df_merge$AnnotatedID, df_merge$TimePoint, sep = ".")
df_selected = df_merge[,c("id", "id_order", "AnnotatedID", "count", "Outcome", "SUBJECT.WUTX")]
h = median(df_selected$count)

jpeg("LgTxCF12_microbial_readCounts.jpg", res = 300, height = 15, width = 30, units = 'cm')
ggplot(df_selected, aes(x = id_order, y = count, fill = AnnotatedID)) + 
  geom_bar(stat="identity", position="stack") +
  ggtitle("# microbial reads per sample") +  
  scale_fill_manual(values = cbbPalette_12) +
  #scale_fill_brewer(palette = "Set3")
  labs(y = "# reads", x = "Sample") + theme_bw() +
  #scale_y_continuous(breaks = c(0, 1)) + 
  theme(axis.text.x = element_text(colour="black", size=10, angle=90, hjust=1, vjust=0.5, face="plain"),
        axis.text.y = element_text(colour="black", size=10, angle=0, hjust=1, vjust=1, face="bold"),
        axis.title.x = element_text(colour="black", size=15, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=15, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=15, face="plain"), 
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "bottom",
        strip.text = element_text(size=12, face = "bold")) +
  guides(colour = guide_legend(nrow = 2)) +
  geom_hline(yintercept=median(df_selected$count)) +
  geom_text(aes(0, h, label = h, vjust = 0, hjust = 0))
dev.off()  



############################################
########## Beta Diversity     ##############
############################################
### TODO: Draw heatmap


############################################
########## Ordination Methods ##############
############################################
## TODO: Visualize 3D graphs
## Calcualte pairwize distances for each coordination
ord_nmds_bray = ordinate(microbial, method = "NMDS", distance = "bray")
ord_nmds_jaccard = ordinate(microbial, method = "NMDS", distance = "jaccard")
ord_nmds_jsd = ordinate(microbial, method = "NMDS", distance = "jsd")

ord_pcoa_bray = ordinate(microbial, method = "PCoA", distance = "bray")
ord_pcoa_jaccard = ordinate(microbial, method = "PCoA", distance = "jaccard")
ord_pcoa_jsd = ordinate(microbial, method = "PCoA", distance = "jsd")


###############################
### Scree Plot: show the proportion of variance that is represented for each axis
################################
plot_scree(ord_pcoa_bray, "Scree Plot: Bray-Curtis PCoA")
plot_scree(ord_pcoa_jaccard, "Scree Plot: Jaccard PCoA")

jpeg("LgTxCF12_screePlot_pcoa_jsd.jpg", res = 300, height = 10, width = 20, units = 'cm')
plot_scree(ord_pcoa_jsd) +
ggtitle("Scree Plot: JSD PCoA") +  
  scale_fill_manual(values = cbbPalette_12) +
  scale_color_manual(values=c( "red", "midnightblue"),  breaks=c("bos", "non-bos")) +
  #scale_fill_brewer(palette = "Set3")
  #labs(y = "# reads", x = "Sample") + theme_bw() +
  #scale_y_continuous(breaks = c(0, 1)) + 
  theme_bw() + 
  theme(axis.text.x = element_text(colour="black", size=7, angle=90, hjust=1, vjust=1, face="plain"),
        axis.text.y = element_text(colour="black", size=10, angle=0, hjust=1, vjust=0.5, face="plain"),
        axis.title.x = element_text(colour="black", size=15, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=15, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=15, face="plain"), 
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "top",
        strip.text = element_text(size=12, face = "bold")) +
  guides(colour = guide_legend(nrow = 1))
dev.off()

######################################################
## NMDS (non-parametric Multi-Dimentional Scaling)  ##
######################################################
jpeg("LgTxCF12_ordination_nmds_jaccard_v0.2.jpg", res = 300, height = 15, width = 15, units = 'cm')
plot_ordination(microbial, ord_nmds_jaccard, color = "Outcome") +  geom_point(size = 3) +
  ggtitle("NMDS using Jaccard") +  
  scale_fill_manual(values = cbbPalette_12) +
  scale_color_manual(values=c( "red", "midnightblue"),  breaks=c("BO", "non-BO")) +
  #scale_fill_brewer(palette = "Set3")
  #labs(y = "# reads", x = "Sample") + theme_bw() +
  #scale_y_continuous(breaks = c(0, 1)) + 
  theme_bw() + 
  theme(axis.text.x = element_text(colour="black", size=10, angle=0, hjust=0.5, vjust=0.5, face="plain"),
        axis.text.y = element_text(colour="black", size=10, angle=0, hjust=1, vjust=0.5, face="plain"),
        axis.title.x = element_text(colour="black", size=15, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=15, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=15, face="plain"), 
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "top",
        strip.text = element_text(size=12, face = "bold")) +
  guides(colour = guide_legend(nrow = 1))
dev.off() 




jpeg("LgTxCF12_ordination_nmds_jaccard_Intapatient_v0.1.jpg", res = 300, height = 15, width = 25, units = 'cm')
plot_ordination(microbial, ord_nmds_jaccard, color = "SUBJECT.WUTX") +  geom_point(size = 3) +
  ggtitle("NMDS using Jaccard") +  
  scale_fill_manual(values = cbbPalette_12) +
  #scale_color_manual(values=c( "red", "midnightblue"),  breaks=c("BO", "non-BO")) +
  #scale_fill_brewer(palette = "Set3")
  #labs(y = "# reads", x = "Sample") + theme_bw() +
  #scale_y_continuous(breaks = c(0, 1)) + 
  theme_bw() + 
  theme(axis.text.x = element_text(colour="black", size=10, angle=0, hjust=0.5, vjust=0.5, face="plain"),
        axis.text.y = element_text(colour="black", size=10, angle=0, hjust=1, vjust=0.5, face="plain"),
        axis.title.x = element_text(colour="black", size=15, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=15, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=15, face="plain"), 
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "top",
        strip.text = element_text(size=12, face = "bold")) +
  guides(colour = guide_legend(nrow = 1))
dev.off() 




jpeg("LgTxCF12_ordination_nmds_jaccard_timepoint_v0.1.jpg", res = 300, height = 15, width = 25, units = 'cm')
plot_ordination(microbial, ord_nmds_jaccard, color = "TimeAfterTransplantation") +  geom_point(size = 3) +
  ggtitle("NMDS using Jaccard") +  
  scale_fill_manual(values = cbbPalette_12) +
  #scale_color_manual(values=c( "red", "midnightblue"),  breaks=c("BO", "non-BO")) +
  #scale_fill_brewer(palette = "Set3")
  #labs(y = "# reads", x = "Sample") + theme_bw() +
  #scale_y_continuous(breaks = c(0, 1)) + 
  theme_bw() + 
  theme(axis.text.x = element_text(colour="black", size=10, angle=0, hjust=0.5, vjust=0.5, face="plain"),
        axis.text.y = element_text(colour="black", size=10, angle=0, hjust=1, vjust=0.5, face="plain"),
        axis.title.x = element_text(colour="black", size=15, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=15, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=15, face="plain"), 
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "top",
        strip.text = element_text(size=12, face = "bold")) +
  guides(colour = guide_legend(nrow = 1))
dev.off() 




###################
## PCoA = MDS  ####
###################
jpeg("LgTxCF12_ordination_pcoa_jsd.jpg", res = 300, height = 15, width = 15, units = 'cm')
plot_ordination(microbial, ord_pcoa_jsd, color = "SUBJECT.WUTX") +  geom_point(size = 3) +
  ggtitle("PCoA using JSD") +  
  scale_fill_manual(values = cbbPalette_12) +
  #scale_color_manual(values=c( "red", "midnightblue"),  breaks=c("bos", "non-bos")) +
  #scale_fill_brewer(palette = "Set3")
  #labs(y = "# reads", x = "Sample") + theme_bw() +
  #scale_y_continuous(breaks = c(0, 1)) + 
  theme_bw() + 
  theme(axis.text.x = element_text(colour="black", size=10, angle=0, hjust=0.5, vjust=0.5, face="plain"),
        axis.text.y = element_text(colour="black", size=10, angle=0, hjust=1, vjust=0.5, face="plain"),
        axis.title.x = element_text(colour="black", size=15, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=15, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=15, face="plain"), 
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "top",
        strip.text = element_text(size=12, face = "bold")) +
  guides(colour = guide_legend(nrow = 1))
dev.off() 




#################
### t-SNE  ######
#################
mic_uniq = t(unique(otu_table(microbial)))
tsne_microbial <- Rtsne(otu_table(mic_uniq), dim = 2, perplexity=20, 
                        verbose=TRUE, max_iter = 2000)

## TODO: sample S11.E is not ordered, it is located on the last row 
tsne_df = data.frame(tsne = tsne_microbial$Y, outcome = sample_data(microbial)$Outcome, subject = sample_data(microbial)$SUBJECT.WUT)
jpeg("LgTxCF12_tsne.jpg", res = 300, height = 15, width = 15, units = 'cm')
ggplot(tsne_df, aes(x = tsne.1, y= tsne.2, color = subject)) +
  geom_point(size = 3) +
  ggtitle("t-SNE") +  
  scale_fill_manual(values = cbbPalette_12) +
  scale_color_manual(values = cbbPalette_12) +
#  scale_color_manual(values=c( "red", "midnightblue"),  breaks=c("bos", "non-bos")) +
  #scale_fill_brewer(palette = "Set3")
  #labs(y = "# reads", x = "Sample") + theme_bw() +
  #scale_y_continuous(breaks = c(0, 1)) + 
  theme_bw() + 
  theme(axis.text.x = element_text(colour="black", size=10, angle=0, hjust=0.5, vjust=0.5, face="plain"),
        axis.text.y = element_text(colour="black", size=10, angle=0, hjust=1, vjust=0.5, face="plain"),
        axis.title.x = element_text(colour="black", size=15, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=15, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=15, face="plain"), 
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "top",
        strip.text = element_text(size=12, face = "bold")) +
  guides(colour = guide_legend(nrow = 2))
dev.off() 

