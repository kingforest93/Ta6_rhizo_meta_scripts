#clean
rm(list=ls())
gc()
setwd("D:/项目试验方案数据/小麦微生物组/指标数据/宏基因组分析结果/06.function")

#load package
library(psych)
library(pheatmap)
library(RColorBrewer)
library(ggtree)

#read data
tax <- read.csv("combined_42.species_abundance", header=T, sep="\t", row.names=1)
tax <- as.data.frame(t(tax))
soil <- read.csv("soil_metadata.csv", header=T, sep=",")
plant <- read.csv("plant_metadata.csv", header=T, sep=",")
#HZn_vs_LZn <- read.csv("OG_wilcox_YL_HZn_vs_LZn.csv", header=T)
#HZn_vs_Soil <- read.csv("OG_wilcox_YL_HZn_vs_Soil.csv", header=T)
#LZn_vs_Soil <- read.csv("OG_wilcox_YL_LZn_vs_Soil.csv", header=T)
HLZn <- read.csv("species_wilcox_YL_HZn_vs_LZn.csv", header=T)

#correlation between taxa and metadata
Rhizo <- read.csv("species_wilcox_YL_Rhizo_vs_Soil_nomin.csv", header=T)
#KO_ZSB <- read.csv("KEGG_KO_Zn_mobilization.txt", header=F)
#KO_ZSB <- as.vector(KO_ZSB[,1])
KO_ZSB <- intersect(Rhizo[Rhizo$Group=='Rhizo','Taxon'], 
                    HLZn[HLZn$Group!="", 'Taxon'])
#Rhizo <- intersect(HZn_vs_Soil[HZn_vs_Soil$Group=="HZn",'Taxon'],
#                   LZn_vs_Soil[LZn_vs_Soil$Group=="LZn",'Taxon'])
#Diff <- c(intersect(Rhizo, HZn_vs_LZn[HZn_vs_LZn$Group=="HZn",'Taxon']),
#          intersect(Rhizo, HZn_vs_LZn[HZn_vs_LZn$Group=="LZn",'Taxon']))
#tax.dat <- tax[grepl("YS.[1-4].[^B]", rownames(tax)), colnames(tax) %in% Diff]
#tax.dat <- tax[grepl("YS.[1-4].[^B]", rownames(tax)), colnames(tax) %in% Rhizo]
tax.dat <- tax[grepl("YL.[1-4].[^B]", rownames(tax)), colnames(tax) %in% KO_ZSB]
tax.dat <- tax.dat[,order(colnames(tax.dat))]
#tax.dat <- tax.dat[,order(apply(tax.dat, 2, function(x)max(x, na.rm=TRUE)),
#                          decreasing=TRUE)]
soil.dat <- soil[order(soil$Block),]
soil.dat <- subset(soil.dat, Site=="YL" & Cultivar!="bulk", c('DTPA_Zn'))
plant.dat <- plant[order(plant$Block),]
plant.dat <- subset(plant.dat, Site=="YL", c('AnPlZnU', 'MtPlZnU', 'MtGrZnC'))
metadata <- cbind(soil.dat, plant.dat)
rownames(metadata) <- rownames(tax.dat)
cor.mx <- corr.test(x=tax.dat, y=metadata, method="spearman", adjust="BH", ci=FALSE)
#cor.mx$p <- cor.mx$p[apply(cor.mx$r, 1, max)>0.5 | apply(cor.mx$r, 1, min)<(-0.5),]
#cor.mx$r <- cor.mx$r[apply(cor.mx$r, 1, max)>0.5 | apply(cor.mx$r, 1, min)<(-0.5),]
write.csv(cor.mx$r, "species_rhizo_HLZn_metadata_correlation_YL.csv", row.names=TRUE)
pheatmap(cor.mx$r, color=colorRampPalette(c(brewer.pal(3, "Set2")[1], "white", 
                                            brewer.pal(3, "Set2")[2]))(100),
         display_numbers=ifelse(cor.mx$p<0.05,'*',''), cellwidth=15,
         cellheight=10, main="Spearman correlation", annotation_row=NULL, 
         annotation_colors=NULL, annotation_col=NULL,
         cluster_rows=FALSE, cluster_cols=FALSE,
         filename="species_rhizo_HLZn_metadata_corrplot_YL.pdf",
         width=7, height=44)

#dendrogram of selected taxa
lineage <- read.table("genus.lineage", header=F)
cor.mx$r <- cor.mx$r[apply(cor.mx$r, 1, min)<(-0.5) | apply(cor.mx$r, 1, max)>0.5,]
tax.sel <- data.frame(Taxon=rownames(cor.mx$r), Abund=rep(0.0, nrow(cor.mx$r)),
                      Lineage=rep('', nrow(cor.mx$r)))
for (i in 1:nrow(cor.mx$r)) {
  tax.sel$Lineage[i] <- lineage[grep(tax.sel$Taxon[i], lineage$V1),]
  tax.sel$Abund[i] <- median(tax.dat[,tax.sel$Taxon[i]])
}

#small test
t <- read.tree(text="((((((S1,S2)Genus)Family)Order)Class)Phylum)Root;")
ggtree(tr=t, layout="rectangular", branch.length="none") +
  geom_tiplab() +
  geom_tippoint() +
  geom_nodepoint() +
  geom_nodelab()






