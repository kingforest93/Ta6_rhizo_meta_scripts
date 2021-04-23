#clean
rm(list=ls())
gc()
setwd(dir="D:/10-AGIS博后/项目试验方案数据/小麦微生物组/指标数据/宏基因组分析结果/06.function")

#load package
library(vegan)
library(agricolae)
library(tidyr)
library(reshape2)
library(ggplot2)
library(ggsignif)
library(pheatmap)
library(RColorBrewer)

#read tax count matrix
tax <- read.csv("combined_42.genus_abundance", sep="\t", header=TRUE, row.names=1)
rhizo <- read.csv("genus_wilcox_YL_Rhizo_vs_Soil.csv", header=TRUE)
rhizo <- rhizo[rhizo$Group=="Rhizo", 'Taxon']
tax_dat <- tax[rownames(tax) %in% rhizo, grepl("YL.[1-4].[^B]", colnames(tax))]
# tax_dat <- tax[, grepl("YL", colnames(tax))]
tax_dat <- as.matrix(t(tax_dat))

#α-diversity boxplot
tax_div <- diversity(tax_dat, index="shannon")
tax_div <- data.frame(sample=names(tax_div), shannon=as.vector(tax_div))
# tax_div$group <- as.factor(rep(c(rep('Rhizo', 6), 'Soil'), 4))
tax_div$group <- as.factor(rep(c('LZn','HZn','LZn','LZn','HZn','HZn'),4))
# tax_div$group <- as.factor(substr(rownames(tax_dat), 6, 7))
tax_test <- t.test(shannon~group, tax_div, var.equal=FALSE, paired=FALSE)
sig <- ifelse(tax_test$p.value<0.001, '***', ifelse(tax_test$p.value<0.01, '**',
              ifelse(tax_test$p.value<0.05, '*', 'ns.')))
tax_aov <- data.frame(group=levels(tax_div$group),
                      mean=tapply(tax_div$shannon, tax_div$group, mean),
                      se=tapply(tax_div$shannon, tax_div$group, function(x)
                        sqrt(sd(x) / length(x))))
top <- max(tax_aov$mean)+max(tax_aov$se)
#tax_aov <- aov(shannon~group, tax_div)
#div_aov <- duncan.test(tax_aov, 'group', alpha=0.05, group=TRUE)
#div_aov$groups <- div_aov$groups[c('3','4','14','10','17','18','B'),]
#div_aov$means <- div_aov$means[c('3','4','14','10','17','18','B'),]
#tax_aov <- data.frame(group=factor(rownames(div_aov$means),
#                                   levels=c('3','4','14','10','17','18','B'), 
#                                   ordered=TRUE),
#                      mean=div_aov$means$shannon,
#                      se=div_aov$means$std / sqrt(div_aov$means$r),
#                      sig=div_aov$groups$groups)
bar.div <- ggplot(tax_aov, aes(x=group, y=mean, ymin=mean-se, ymax=mean+se)) +
  geom_col(aes(fill=group), color="black", width=0.6, size=0.3) +
  geom_errorbar(color="black", width=0.15, size=0.3) +
  annotate("text", x=1.5, y=top+0.4, label=sig, size=3) +
  annotate("segment", x=1, xend=2, y=top+0.2, yend=top+0.2, color="black") +
  #geom_text(aes(x=group, y=mean+se+0.2, label=sig), size=3) +
  labs(x=NULL, y="OG shannon index", fill=NULL) +
  scale_y_continuous(breaks=seq(0, 5, 1), limits=c(0, 5)) +
  scale_fill_manual(values=brewer.pal(12, "Paired")[c(6,2)]) +
  # scale_fill_manual(values=brewer.pal(12, "Paired")[c(4,8)]) +
  theme(line=element_line(color="black", size=0.3),
        text=element_text(color="black", size=8),
        axis.ticks=element_line(color="black", size=0.3),
        axis.text=element_text(color="black", size=8),
        panel.background=element_rect(fill="white", color="black"),
        panel.grid=element_blank(), legend.position="none")
box.div <- ggplot(tax_div, aes(x=group, y=shannon)) +
  geom_boxplot(aes(fill=group), width=0.6, color="black", 
               size=0.3, outlier.size=0.5) +
  geom_signif(comparisons=list(c('Rhizo', 'Soil')), #comparisons=list(c('Rhizo', 'Soil')),
              test="t.test", map_signif_level=TRUE, tip_length=0.02, 
              step_increase=0.05, textsize=2, size=0.3) +
  labs(x=NULL, y="KO shannon index", fill=NULL) +
  #scale_fill_manual(values=brewer.pal(12, "Paired")[c(6,2)]) +
  scale_fill_manual(values=brewer.pal(12, "Paired")[c(4,8)]) +
  theme(line=element_line(color="black", size=0.3), 
        text=element_text(color="black", size=8),
        axis.ticks=element_line(color="black", size=0.3),
        axis.text=element_text(color="black", size=8),
        panel.background=element_rect(fill="white", color="black"),
        panel.grid=element_blank(), legend.position="none")
ggsave("rhizo_genus_shannon_barplot_HZn_LZn_YL.pdf", bar.div, width=4, height=6, units="cm")

#β-diversity heatmap and clustering
tax.dist <- vegdist(tax_dat, method="bray", diag=TRUE, upper=TRUE)
tax.dist <- as.matrix(tax.dist)
annot.row <- data.frame(Group=as.factor(rep(c('LZn', 'HZn', 'LZn', 'LZn', 'HZn', 'HZn', 'Soil'), 4)),
                        row.names=rownames(tax.dist))
cols <- brewer.pal(12, "Paired") 
annot.color <- list(Group=c(HZn=cols[6], LZn=cols[2], Soil=cols[4]))
pheatmap(tax.dist, clustering_method="ward.D", 
         color=colorRampPalette(c("grey", "white"))(100),
         display_numbers=TRUE, cellwidth=30, cellheight=18,
         main="Bray-Curtis dissimilarity", annotation_row=annot.row, 
         annotation_colors=annot.color, annotation_col=NULL,
         legend=FALSE, show_colnames=FALSE, treeheight_col=0, treeheight_row=100,
         cutree_rows=3, cutree_cols=3, filename="OG_bray_heamap_YL.pdf")

#stacked barplot
tax <- read.csv("combined_42.phylum_abundance", sep="\t", header=TRUE, row.names=1)
#rhizo <- read.csv("phylum_wilcox_YL_Rhizo_vs_Soil_nomin.csv", header=TRUE)
#rhizo <- rhizo[rhizo$Group=="Rhizo", 'Taxon']
#tax <- tax[rownames(tax) %in% rhizo, grepl("YL.[1-4].[^B]", colnames(tax))]
tax <- tax[, grepl("YL", colnames(tax))]
tax <- as.data.frame(t(tax))
tax_major <- tax[,colMeans(tax)>=1]
tax_minor <- tax[,colMeans(tax)<1]
tax_major$ZOther <- rowSums(tax_minor)
group <- as.factor(rep(c(rep('Rhizo', 6),'Soil'), 4))
#group <- as.factor(rep(c('LZn','HZn','LZn','LZn','HZn','HZn'),4))
#group <- factor(substr(rownames(tax_major), 6, 7))
tmp <- tax_major[1:2,]
#tmp <- tax_major[1:7,]
rownames(tmp) <- levels(group)
for (i in 1:ncol(tax_major)) {
  tmp[,i] <- tapply(tax_major[,i], group, mean)
}
tmp$group <- factor(rownames(tmp), 
                    levels=c('Rhizo', 'Soil'),
                    #levels=c('HZn','LZn'),
                    #levels=c('3','4','14','10','17','18','B'),
                    ordered=TRUE)
tax_major <- tmp
tax_major <- gather(tax_major, key="phylum", value="abundance", -group)
tax.p <- ggplot(tax_major, aes(x=group, y=abundance)) +
  geom_bar(aes(fill=phylum), stat="identity", position="stack", 
           color="black", width=0.6, size=0.3) +
  scale_fill_manual(values=brewer.pal(12, "Paired")) +
  labs(x=NULL, y="Relative abundance (%)", fill="Phylum") +
  theme(line=element_line(color="black", size=0.3), 
        text=element_text(color="black", size=8),
        axis.ticks=element_line(color="black", size=0.3),
        axis.text=element_text(color="black", size=8),
        panel.background=element_rect(fill="white", color="black"),
        panel.grid=element_blank(), legend.key=element_rect(fill="white"))
ggsave("phylum_stacked_barplot_Rhizo_Soil_YL.pdf", tax.p, width=8, height=8, units="cm")

#boxplot with wilcox test
tax <- read.csv("combined_42.COG_GROUP_abundance", sep="\t", header=TRUE, row.names=1)
rhizo <- read.csv("COG_wilcox_YL_Rhizo_vs_Soil_nomin.csv", header=TRUE)
rhizo <- rhizo[rhizo$Group=="Rhizo", 'Taxon']
tax <- tax[rownames(tax) %in% rhizo, grepl("YL.[1-4].[^B]", colnames(tax))]
tax <- tax[order(apply(tax, 1, function(x)max(x, na.rm=TRUE)), decreasing=TRUE),]
tax <- as.data.frame(t(tax))
ord <- colnames(tax)
#tax$group <- as.factor(rep(c(rep('Rhizo', 6), 'Soil'), 4))
tax$group <- as.factor(rep(c('LZn','HZn','LZn','LZn','HZn','HZn'),4))
anno <- data.frame(phylum=rep(ord, each=3), group=rep(c("HZn","LZn","Soil"), length(ord)),
                  pos=rep(rep(0,3), length(ord)), max=rep(rep(0,3), length(ord)), 
                  sig=rep(rep("",3), length(ord)))
for (i in 1:length(ord)) {
  tmp <- kruskal(tax[,i], tax$group, p.adj="BH", alpha=0.05)
  anno$max[(3*(i-1)+1):(3*i)] <- tmp$means[c('HZn','LZn','Soil'),'Max']
  anno$sig[(3*(i-1)+1):(3*i)] <- tmp$groups[c('HZn','LZn','Soil'),'groups']
  anno$pos[(3*(i-1)+1):(3*i)] <- c(i-0.2,i,i+0.2)
}
#tax$group <- factor(substr(rownames(tax),6,7), 
#                    levels=c('3','4','14','10','17','18','B'), ordered=TRUE)
#anno <- data.frame(phylum=rep(ord, each=7), 
#                   group=rep(c('3','4','14','10','17','18','B'), length(ord)),
#                   pos=rep(rep(0,7), length(ord)), max=rep(rep(0,7), length(ord)), 
#                   sig=rep(rep("",7), length(ord)))
#for (i in 1:length(ord)) {
#  tmp <- kruskal(tax[,i], tax$group, p.adj="BH", alpha=0.05)
#  anno$max[(7*(i-1)+1):(7*i)] <- tmp$means[c('3','4','14','10','17','18','B'),'Max']
#  anno$sig[(7*(i-1)+1):(7*i)] <- tmp$groups[c('3','4','14','10','17','18','B'),'groups']
#  anno$pos[(7*(i-1)+1):(7*i)] <- c('3','4','14','10','17','18','B')
#}
tax <- gather(tax, key="phylum", value="abundance", -group)
tax$phylum <- factor(x=tax$phylum, levels=unique(tax$phylum), ordered=TRUE)
#anno$phylum <- factor(x=anno$phylum, levels=unique(anno$phylum), ordered=TRUE)
#tax.p <- ggplot(tax, aes(x=phylum, y=abundance)) +
#  geom_boxplot(aes(fill=group), color="black", width=0.6, size=0.3,
#               outlier.size=0.3) +
#  geom_text(data=anno, aes(x=pos, y=max+0.1, label=sig), size=1.5) +
#  scale_fill_manual(values=brewer.pal(12, "Paired")[c(6,6,6,2,2,2,4)]) +
#  labs(x="COG_GROUP", y="Relative abundance (%)", fill=NULL) +
#  theme(line=element_line(color="black", size=0.3), 
#        text=element_text(color="black", size=8),
#        axis.ticks=element_line(color="black", size=0.3),
#        axis.text=element_text(color="black", size=8),
#        axis.text.x=element_text(color="black", size=8, angle=45, vjust=0.5),
#        panel.background=element_rect(fill="white", color="black"),
#        panel.grid=element_blank(), legend.key=element_rect(fill="white"),
#        legend.position="top")
tax.p <- ggplot(tax, aes(x=group, y=abundance)) +
  geom_boxplot(aes(fill=group), color="black", width=0.6, size=0.3, 
               outlier.size=0.3) +
#  geom_text(data=anno, aes(x=pos, y=max*1.05, label=sig), size=1.5) +
  geom_signif(comparisons=list(c('HZn', 'LZn')),
              #comparisons=list(c('Rhizo', 'Soil')),
              test="wilcox.test", map_signif_level=TRUE, tip_length=0.02, 
              step_increase=0.05, textsize=2, size=0.3) +
  scale_fill_manual(values=brewer.pal(12, "Paired")[c(6, 2)]) +
  labs(x="COG_GROUP", y="Relative abundance (%)", fill=NULL) +
  facet_wrap(~phylum, nrow=1, strip.position="bottom", scales="fixed") +
  theme(line=element_line(color="black", size=0.3), 
        text=element_text(color="black", size=8),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        panel.background=element_rect(fill="white", color="black"),
        panel.grid=element_blank(), legend.key=element_rect(fill="white"),
        legend.position="top")
ggsave("rhizo_COG_boxplot_HZn_LZn_YL.pdf", tax.p, width=5, height=6, units="cm")
