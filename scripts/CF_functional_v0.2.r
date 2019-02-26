require("heatmap3")
library("data.table")
library(ggplot2)
library("phyloseq")
library("DESeq2")
require(vegan)
#library("Rtsne")
dev.off()
rm(list=ls())



setwd("/Users/ahmedmetwally/Dropbox/LungTransplantWASHU")
setwd("D:/Dropbox/LungTransplantWASHU")

# The palette with grey:
cbPalette  = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbbPalette = c("#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000", "#E69F00", "#56B4E9", "#009E73")
cbPalette  = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbbPalette_10 = c("#999999", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000", "#E69F00", "#56B4E9", "#009E73", "darkorchid")
cbbPalette_9 = c("#56B4E9", "sienna1", "#009E73", "darkorchid",  "#F0E442", "#0072B2", "#CC79A7", "#999999", "#000000")
cbbPalette_12 = c("#999999", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000", "#E69F00", "#56B4E9", "#009E73", "darkorchid",
                  "mistyrose1", "sienna1")



### Prepare count matrix. Filter taxa that have read count < 5 reads 
funcTable = read.csv(file="humann2/WASHU_CF_LG_TX_genefamilies_cpm.csv", header=TRUE, check.names = FALSE, row.names = 1)
tmp = as.matrix(funcTable)



# row.names(tmp)=tmp[,1]
# tmp=tmp[,-1]
funcTable_2 = apply(tmp, 2, function(x){
  #x[(x/sum(x)) < 0.01] = 0 
  x[x < 5] = 0 
  x
})


funcTable_3 = funcTable_2[rowSums(funcTable_2)!=0,]
funcTable_3 = funcTable_3[!grepl("|", rownames(funcTable_3), fixed=TRUE),] ## Remove stratified info



## Prepare Taxa matrix
# taxaTable = read.csv(file="taxProfiles/CF12/LungTxCF12_OTUMatrix.csv", header=TRUE, check.names = FALSE, row.names = 1)
# taxaTable_2 = as.matrix(taxaTable)

# Prepare metadata matrix
meta = read.csv(file="metadata/CF_12_v1.7.csv", header = TRUE)



## Draw Histogram
sort(apply(funcTable, 1, sum),  decreasing = TRUE)[1:20]
hist(as.integer(funcTable["UniRef90_O00370: LINE-1 retrotransposable element ORF2 protein|unclassified",]))


## Man-Whittney
meta_sorted = meta[with(meta, order(meta$SampleIDFP)), ]
func_mw = t(apply(funcTable_3, 1, function(x) unlist(wilcox.test(x ~ meta_sorted$Outcome))))
func_mw_df = as.data.frame(func_mw)
func_mw_sorted = func_mw_df[with(func_mw_df, order(p.value)), ]
write.csv(func_mw_sorted, "Functional_mannwhittney.csv")
top10sig = funcTable_3[rownames((func_mw_sorted))[1:10],]
rownames(top10sig) = gsub(":", "_", rownames(top10sig), fixed=TRUE)
rownames(top10sig) = gsub(" ", "_", rownames(top10sig), fixed=TRUE)
rownames(top10sig) = gsub("-", "_", rownames(top10sig), fixed=TRUE)
rownames(top10sig) = gsub("(", "_", rownames(top10sig), fixed=TRUE)
rownames(top10sig) = gsub(")", "_", rownames(top10sig), fixed=TRUE)



## t-test

func_ttest = t(apply(funcTable_3, 1, function(x) unlist(t.test(x ~ meta_sorted$Outcome))))
func_ttest_df = as.data.frame(func_ttest)
func_ttest_sorted = func_ttest_df[with(func_ttest_df, order(p.value)), ]
write.csv(func_ttest_sorted, "Functional_ttest.csv")
top10sig = funcTable_3[rownames((func_ttest_sorted))[1:10],]

# MetaLonDA
group.vec = meta_sorted$Outcome
id.vec = meta_sorted$SUBJECT.WUTX
time.vec = meta_sorted$TimeAfterTransplantation


count.matrix = top10sig

library(MetaLonDA)
output_all_nbinomial_2 = metalondaAll(Count = count.matrix, Time = time.vec, Group = group.vec, 
                                      ID = id.vec, fit.method = "lowess", n.perm = 1000, 
                                      num.intervals = 99, parall = FALSE, pvalue.threshold = 0.05, 
                                      adjust.method = "none", time.unit = "days", norm.method = "none",
                                      prefix = "Lx", col = c("firebrick", "blue"), ylabel = "Normalized count (cpm)")





##### Visualization



#meta_trans = t(meta)
top10sig_trans = t(top10sig)
top10sig_trans = as.data.frame(top10sig_trans)
top10sig_trans$SampleIDFP=rownames(top10sig_trans)
meta_top10sig = merge(meta,top10sig_trans)




x="UniRef90_A9AHG3__Sulfate_adenylyltransferase_subunit_2"
jpeg(paste("LgTxCF12_functional_", x, ".jpg", sep=""), res = 300, height = 20, width = 30, units = 'cm')
p = ggplot(meta_top10sig, aes(x = SampleIDFP, y = UniRef90_A9AHG3__Sulfate_adenylyltransferase_subunit_2, fill = Outcome)) + theme_bw() 
p + geom_bar(stat="identity", position="stack") + theme(legend.position="bottom") + 
  ggtitle(x) + 
  labs(y = "Normalized count (cpm)", x = "Sample") +
  #scale_fill_manual(values = cbbPalette_12) +
  scale_fill_manual(values=c( "red", "midnightblue"),  breaks=c("bos", "non-bos")) + 
  scale_color_manual(values=c( "red", "midnightblue"),  breaks=c("bos", "non-bos")) +
  theme(axis.text.x = element_text(colour="black", size=10, angle=90, hjust=0.5, vjust=0, face="bold"),
        axis.text.y = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="bold"),
        axis.title.x = element_text(colour="black", size=15, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=15, angle=90, hjust=.5, vjust=.5, face="bold"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.text=element_text(size=15, face="plain"), legend.title = element_blank())
dev.off()


colnames(top10sig_trans)


