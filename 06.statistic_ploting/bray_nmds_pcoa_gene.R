rm(list=ls())
gc()
setwd(dir="D:/10-AGIS博后/项目试验方案数据/小麦微生物组/指标数据/宏基因组分析结果/06.function")

#load package
library(vegan)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)

#read taxonomy count matrix
tax <- read.csv("combined_42.species_abundance", sep="\t", header=TRUE, row.names=1)
tax_dat <- tax[,grepl("YL", colnames(tax))]
# rhizo <- read.csv("species_wilcox_YL_Rhizo_vs_Soil.csv", header=TRUE)
# rhizo <- rhizo[rhizo$Group=="Rhizo", 'Taxon']
# tax_dat <- tax[rownames(tax) %in% rhizo, grepl("YL.[1-4].[^B]", colnames(tax))]
tax_dat <- data.frame(t(tax_dat))

#PCA analysis
#tax.pca <- rda(tax_dat, scale=FALSE)
#var.perc <- round(tax.pca$CA$eig[1:2] / tax.pca$CA$tot.chi * 100, 1)
#plot(tax.pca, display="si", choices=c(1,2), 
#     main="PCA by SVD using vegan-rda",
#     xlab=paste("PC1 (", var.perc[1], "%)", sep=""),
#     ylab=paste("PC1 (", var.perc[2], "%)", sep=""))
#pca.score <- as.data.frame(tax.pca$CA$u[,1:2])

#PCoA analysis
tax.dist <- vegdist(tax_dat, method="bray", diag=TRUE, upper=TRUE)
tax_pcoa <- cmdscale(tax.dist, eig=TRUE)
var.perc <- tax_pcoa$eig[1:2] / sum(tax_pcoa$eig) * 100
#ordiplot(tax_pcoa, type = "text", display="sites", main="PCoA/Bray", 
#	xlab=paste("PCoA1 (", round(var.perc[1],1), "%)", sep=""), 
#	ylab=paste("PCoA2 (", round(var.perc[2],1), "%)", sep=""))

#draw PCoA plot with ggplot and ggrepel
#tax.clus <- as.factor(rep(c('LZn','HZn','LZn','LZn','HZn','HZn'),4))
tax.clus <- as.factor(rep(c(rep('Rhizo', 6), 'Soil'), 4))
#tax_adonis <- adonis(tax.dist~tax.clus)
#adonis.F <- tax_adonis$aov.tab$F.Model[1]
#adonis.sig <- tax_adonis$aov.tab$'Pr(>F)'[1]
pcoa.score <- as.data.frame(scores(tax_pcoa))
colnames(pcoa.score) <- c('PCoA1', 'PCoA2')
#colnames(pca.score) <- c('PCA1', 'PCA2')
pcoa.p <- ggplot(data=pcoa.score, aes(x=PCoA1, y=PCoA2)) +
  geom_point(aes(color=tax.clus), shape=4, size=1) +
  stat_ellipse(aes(group=tax.clus, color=tax.clus), level=0.9, size=0.2,
               type="norm", linetype=2) +
#  geom_text_repel(aes(label=rownames(pcoa.score)), size=2) +
  labs(x=paste("PCoA1 (", round(var.perc[1],1), "%)", sep=""),
       y=paste("PCoA2 (", round(var.perc[2],1), "%)", sep=""),
       subtitle="Bray-Curtis dissimilarity of species",
       color=NULL) +
#  annotate("text", x=-0.25, y=-0.25, 
#           label=paste("ADONIS\nF = ", round(adonis.F, 3), "\nP = ",
#                       round(adonis.sig, 3), sep=""), size=2.5) +
  #scale_color_manual(values=brewer.pal(12,"Paired")[c(6,2)]) +
  scale_color_manual(values=brewer.pal(12,"Paired")[c(4,8)]) +
  theme(line=element_line(color="black", size=0.3), 
        axis.ticks=element_line(color="black", size=0.3),
        text=element_text(color="black", size=8), 
        axis.text=element_text(color="black", size=8),
        panel.background=element_rect(fill="white", color="black"),
        panel.grid=element_blank(), legend.key=element_rect(fill="white"))
ggsave("species_bray_pcoa_Rhizo_Soil_YL.pdf", pcoa.p, device="pdf", width=10, height=8, units="cm")

#NMDS analysis
tax_dat <- tax_dat * 100
tax_nmds <- metaMDS(tax_dat, distance="bray")
#ordiplot(tax_nmds, type="text", display="sites", main="NMDS/Bray", 
#		sub=paste("Stress:", round(tax_nmds$stress,3)))

#draw NMDS plot with ggplot and ggrepel
#tax.clus <- as.factor(rep(c('LZn','HZn','LZn','LZn','HZn','HZn'),4))
tax.clus <- as.factor(rep(c(rep('Rhizo', 6), 'Soil'), 4))
#tax_anosim <- anosim(tax.dist, tax.clus)
#anosim.R <- tax_anosim$statistic
#anosim.sig <- tax_anosim$signif
nmds.score <- as.data.frame(scores(tax_nmds))
nmds.p <- ggplot(data=nmds.score, aes(x=NMDS1, y=NMDS2)) +
	geom_point(aes(color=tax.clus), shape=4, size=1) +
  #geom_text_repel(aes(label=rownames(nmds.score)), size=2) +
  stat_ellipse(aes(group=tax.clus, color=tax.clus), level=0.9, size=0.2,
               type="norm", linetype=2) +
  labs(x="NMDS1", y="NMDS2",
       subtitle="Bray-Curtis dissimilarity of species", color="Group") +
	#annotate("text", x=-0.27, y=-0.22, 
	#         label=paste("Stress =", round(tax_nmds$stress,3)), size=3) +
  #annotate("text", x=-0.27, y=-0.26, 
  #         label=paste("ANOSIM\nR = ", round(anosim.R, 3), ", P = ", 
  #                     round(anosim.sig, 3), sep=""), size=3) +
  #scale_color_manual(values=brewer.pal(12,"Paired")[c(6,2)]) +
  scale_color_manual(values=brewer.pal(12,"Paired")[c(4,8)]) +
  theme(line=element_line(color="black", size=0.3), 
        axis.ticks=element_line(color="black", size=0.3),
        text=element_text(color="black", size=8), 
        axis.text=element_text(color="black", size=8),
	      panel.background=element_rect(fill="white", color="black"),
	      panel.grid=element_blank(), legend.key=element_rect(fill="white"))
ggsave("species_nmds_bray_Rhizo_Soil_YL.pdf", nmds.p, device="pdf", width=10, height=8, units="cm")
