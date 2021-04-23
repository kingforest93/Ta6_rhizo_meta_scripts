#clean
rm(list=ls())
gc()
setwd("D:/10-AGIS博后/项目试验方案数据/小麦微生物组/指标数据/宏基因组分析结果/06.function")

#load package
library(tidyr)
library(ggplot2)
library(agricolae)
library(ggsignif)
library(ggpubr)
library(RColorBrewer)

#read tax count matrix
tax <- read.csv("combined_42.species_abundance", sep="\t", header=TRUE, row.names=1)
# rhizo <- read.csv("genus_wilcox_YL_Rhizo_vs_Soil.csv", header=TRUE)
# rhizo <- rhizo[rhizo$Group=="Rhizo", 'Taxon']
# tax_dat <- tax[rownames(tax) %in% rhizo, grepl("YL.[1-4].[^B]", colnames(tax))]
#ZSB <- read.csv("OGs_ZSB.txt", header=TRUE, sep="\t")
#tax_dat <- tax[rownames(tax) %in% ZSB$OG, grepl("YL.[1-4].[^B]", colnames(tax))]
tax_dat <- tax[, grepl("YL", colnames(tax))]
tax_dat <- as.data.frame(t(tax_dat))

#wilcoxon test of each taxon
tax_dat$group = factor(rep(c('LZn','HZn','LZn','LZn','HZn','HZn'), 4))
#tax_dat$group = factor(rep(c(rep('Rhizo',6),'Soil'), 4))
dat <- tax_dat[tax_dat$group %in% c('HZn','LZn'),-ncol(tax_dat)]
#dat <- tax_dat[tax_dat$group %in% c('Rhizo','Soil'),-ncol(tax_dat)]
group <- factor(tax_dat[tax_dat$group %in% c('HZn','LZn'),'group'])
#group <- factor(tax_dat[tax_dat$group %in% c('Rhizo','Soil'),'group'])
col.num <- length(colnames(dat))
tax_wilcox <- data.frame(Taxon=colnames(dat), Median=rep(0.0, col.num), 
                         FC=rep(0.0, col.num),
                         P_raw=rep(0.0, col.num), P_adj=rep(0.0, col.num))
tax_wilcox$Median <- apply(dat, 2, median)
tax_wilcox$FC <- apply(dat, 2, function(x) tapply(x, group, median)[1] /
                             tapply(x, group, median)[2])
tax_wilcox$P_raw <- apply(dat, 2, function(x) wilcox.test(x~group, paired=FALSE,
                                                     correct=FALSE)$p.value)
tax_wilcox <- na.omit(tax_wilcox)
tax_wilcox$P_adj <- p.adjust(tax_wilcox$P_raw, method="BH")
tax_wilcox$Group <- rep("", nrow(tax_wilcox))
tax_wilcox[(tax_wilcox$FC>1)&(tax_wilcox$P_adj<0.05), 'Group'] <- "HZn"
#tax_wilcox[(tax_wilcox$FC>1.5)&(tax_wilcox$P_adj<0.05), 'Group'] <- "Rhizo"
tax_wilcox[(tax_wilcox$FC<1)&(tax_wilcox$P_adj<0.05), 'Group'] <- "LZn"
#tax_wilcox[(tax_wilcox$FC<0.5)&(tax_wilcox$P_adj<0.05), 'Group'] <- "Soil"
tax_wilcox$Group <- as.factor(tax_wilcox$Group)
write.csv(tax_wilcox, "rhizo_genus_wilcox_YL_HZn_vs_LZn.csv", row.names=FALSE)
tax_wilcox <- read.csv("ZSB_species_published_potential_median_HZn_vs_LZn.csv", header=TRUE)
tax_wilcox <- tax_wilcox[order(tax_wilcox$Group),]
tax_wilcox.p <- ggplot(tax_wilcox, aes(x=log10(Median), y=log2(FC))) +
                geom_point(aes(color=Group), shape=4, size=1, alpha=1) +
                scale_color_manual(values=c("grey90", brewer.pal(12, "Paired")[c(6,2)])) +
                #scale_color_manual(values=c("grey90", brewer.pal(12, "Paired")[c(4,8)])) +
                labs(x="Log10 (median relative abundance (%))", 
                     y="Log2 (abundance ratio of Rhizo / Soil)",
                     #y="Log2 (abundance ratio of HZn / LZn)",
                     color=NULL) +
                geom_hline(yintercept=log2(c(1.5, 0.5)), linetype="dashed", size=0.3) +
                theme(line=element_line(color="black", size=0.3), 
                      axis.ticks=element_line(color="black", size=0.3), 
                      text=element_text(color="black", size=8),
                      axis.text=element_text(color="black", size=8),
                      panel.background=element_rect(fill="white", color="black"),
                      panel.grid=element_blank(), legend.key=element_rect(fill="white"),
                      legend.text=element_text(color="black", size=8))
x<- ggplot(tax_wilcox, aes(x=log10(Median))) + 
  geom_density(aes(fill=Group), color="black", size=0.3, alpha=0.6) +
  #geom_histogram(aes(fill=Group), bins=100, color="black", size=0.3, alpha=1) +
  scale_fill_manual(values=c("grey90", brewer.pal(12, "Paired")[c(6,2)])) +
  #scale_fill_manual(values=c("grey90", brewer.pal(12, "Paired")[c(4,8)])) +
  theme(line=element_line(color="black", size=0.3), 
        axis.ticks=element_line(color="black", size=0.3), 
        text=element_text(color="black", size=8),
        axis.text=element_text(color="black", size=8),
        panel.background=element_rect(fill="white", color="black"),
        panel.grid=element_blank(), legend.key=element_rect(fill="white"),
        legend.text=element_text(color="black", size=8))
y <- ggplot(tax_wilcox, aes(y=log2(FC))) + 
  geom_density(aes(fill=Group), color="black", size=0.3, alpha=0.6) +
  #geom_histogram(aes(fill=Group), bins=100, color="black", size=0.3, alpha=1) +
  scale_fill_manual(values=c("grey90", brewer.pal(12, "Paired")[c(6,2)])) +
  #scale_fill_manual(values=c("grey90", brewer.pal(12, "Paired")[c(4,8)])) +
  theme(line=element_line(color="black", size=0.3), 
        axis.ticks=element_line(color="black", size=0.3), 
        text=element_text(color="black", size=8),
        axis.text=element_text(color="black", size=8),
        panel.background=element_rect(fill="white", color="black"),
        panel.grid=element_blank(), legend.key=element_rect(fill="white"),
        legend.text=element_text(color="black", size=8))
comb <- ggarrange(x, NULL, tax_wilcox.p, y, nrow=2, ncol=2, 
                  widths=c(3,1), heights=c(1,3), legend="none")
ggsave("rhizo_OG_wilcox_scatter_histo_YL_HZn_vs_LZn.pdf", comb, device="pdf", 
       width=12, height=10, units="cm")
tax_dat <- subset(tax_wilcox, Status=="Potential")
x<- ggplot(tax_dat, aes(x=log10(Median))) + 
  geom_density(aes(fill=Group), color="black", size=0.3, alpha=0.6) +
  #geom_histogram(aes(fill=Group), bins=100, color="black", size=0.3, alpha=1) +
  scale_fill_manual(values=c(brewer.pal(12, "Paired")[c(6,2)])) +
  #scale_fill_manual(values=c(brewer.pal(12, "Paired")[c(4,8)])) +
  labs(x="Log10 (Median relative abundance (%))",
       #x="Log2 (Abundance ratio of Rhizo / Soil)",
       y="Density",
       fill="Group") +
  scale_x_continuous(limits=c(-5, 1), breaks=seq(-5, 1, 1)) +
  theme(line=element_line(color="black", size=0.3), 
        axis.ticks=element_line(color="black", size=0.3), 
        text=element_text(color="black", size=8),
        axis.text=element_text(color="black", size=8),
        panel.background=element_rect(fill="white", color="black"),
        panel.grid=element_blank(), legend.key=element_rect(fill="white"),
        legend.text=element_text(color="black", size=8))
ggsave("ZSB_species_potential_median_density_YL_HZn_vs_LZn.pdf", x, device="pdf",
       width=10, height=5, units="cm")

#kruskal-wallis test of differential taxa
tax <- as.data.frame(t(tax))
HZn_vs_LZn <- read.csv("species_wilcox_YL_HZn_vs_LZn.csv", header=T)
HZn_vs_Soil <- read.csv("species_wilcox_YL_HZn_vs_Soil.csv", header=T)
LZn_vs_Soil <- read.csv("species_wilcox_YL_LZn_vs_Soil.csv", header=T)
Rhizo <- intersect(HZn_vs_Soil[HZn_vs_Soil$Group=="HZn",'Taxon'],
                   LZn_vs_Soil[LZn_vs_Soil$Group=="LZn",'Taxon'])
Diff <- c(intersect(Rhizo, HZn_vs_LZn[HZn_vs_LZn$Group=="HZn",'Taxon']),
          intersect(Rhizo, HZn_vs_LZn[HZn_vs_LZn$Group=="LZn",'Taxon']))
tax_res <- tax[rownames(tax) %in% Diff, grepl("YL", colnames(tax))]
tax_res <- tax_res[order(apply(tax_res, 1, function(x)max(x, na.rm=TRUE)), 
                         decreasing=TRUE),]
tax_res <- as.data.frame(t(tax_res))
tax_res$Group <- as.factor(rep(c('LZn', 'HZn', 'LZn', 'LZn', 'HZn', 'HZn', 'Soil'), 4))
anno_df <- data.frame(Taxon=rep(colnames(tax_res)[-ncol(tax_res)], each=3), 
                      pos=rep(c("HZn", "LZn" , "Soil"), ncol(tax_res)-1),
                      Max=rep(rep(0,3), ncol(tax_res)-1),
                      Sig=rep(rep("",3), ncol(tax_res)-1))
for (i in 1:(ncol(tax_res)-1)) {
  tem <- kruskal(tax_res[,i], tax_res$Group, p.adj="BH", alpha=0.05)
  anno_df$Max[(3*(i-1)+1):(3*i)] <- tem$means[c('HZn','LZn','Soil'),'Max']
  anno_df$Sig[(3*(i-1)+1):(3*i)] <- tem$groups[c('HZn','LZn','Soil'),'groups']
}
anno_df$Taxon <- factor(x=anno_df$Taxon, levels=unique(anno_df$Taxon), ordered=TRUE)
write.csv(anno_df, "species_kruskal_Rhizo_HLZn_YL.csv", row.names=TRUE)
tax_res <- gather(tax_res, key="Taxon", value="Abund", -Group)
tax_res$Taxon <- factor(x=tax_res$Taxon, levels=unique(anno_df$Taxon), ordered=TRUE)
tax_res.p <- ggplot(tax_res, aes(x=Group, y=Abund)) +
  geom_boxplot(aes(fill=Group), color="black", width=0.6, size=0.3, outlier.size=1) +
  geom_text(data=anno_df, aes(x=pos, y=Max+0.04, label=Sig), size=3) +
  facet_wrap(~ Taxon, nrow=17, strip.position="top", scales="fixed") +
  labs(x=NULL, y="Species relative abundance (%)", fill=NULL) +
  scale_fill_manual(values=brewer.pal(12, "Paired")[c(6,2,4)]) +
  theme(line=element_line(color="black", size=0.3), 
        text=element_text(color="black", size=8),
        axis.ticks=element_line(color="black", size=0.3),
        axis.text=element_text(color="black", size=8),
        panel.background=element_rect(fill="white", color="black"),
        panel.grid=element_blank(), legend.key=element_rect(fill="white"),
        legend.position="top")
ggsave("species_kruskal_boxplot_Rhizo_HLZn_YL.pdf", tax_res.p, device="pdf",
       width=34, height=51, units="cm")

#boxplot with wilcox test of ZSB species
tax <- tax_dat
YL_HLZn <- read.csv("rhizo_OG_wilcox_YL_HZn_vs_LZn.csv", header=T)
#YL_HLZn <- YL_HLZn[YL_HLZn$Group!="",]
tax_dat <- tax[, colnames(tax) %in% YL_HLZn$Taxon]
ZSB <- read.csv("OGs_ZSB_selected.txt", header=TRUE, sep="\t")
sel <- intersect(colnames(tax_dat), ZSB$OG)
tax_dat <- as.data.frame(tax_dat[, sel])
rownames(tax_dat) <- rownames(tax)
colnames(tax_dat) <- sel
#tax_dat <- as.data.frame(t(tax_dat))
#tax_dat$Function <- rep("", nrow(tax_dat))
#tax_dat[rownames(tax_dat) %in% YL_ZSB[YL_ZSB$Function=="promoting root Zn uptake",
#                                  'OG'], 'Function'] <- "promoting root Zn uptake"
#tax_dat[rownames(tax_dat) %in% YL_ZSB[YL_ZSB$Function=="Zn ligand metabolism",
#                                  'OG'], 'Function'] <- "Zn ligand metabolism"
# tax_dat[rownames(tax_dat) %in% YL_ZSB[YL_ZSB$Function=="Zn transmembrane transport",
#                                   'OG'], 'Function'] <- "Zn transmembrane transport"
# tax_dat <- apply(tax_dat[,-ncol(tax_dat)], 2, function(x)
#                                 tapply(x, tax_dat$Function, median))
#tax_dat <- as.data.frame(t(tax_dat))
tax_dat <- tax_dat[,order(apply(tax_dat, 2, function(x) max(x, na.rm=TRUE)), 
                         decreasing=TRUE)]
taxon <- factor(x=colnames(tax_dat), 
                levels=unique(colnames(tax_dat)), ordered=TRUE)
# group <- as.factor(rep(c('LZn', 'HZn', 'LZn', 'LZn', 'HZn', 'HZn'), 4))
group <- as.factor(rep(c(rep('Rhizo', 6), 'Soil'), 4))
col.num <- length(colnames(tax_dat))
# tax_wilcox <- data.frame(Taxon=colnames(tax_dat), Median_HZn=rep(0.0, col.num),
#                          Median_LZn=rep(0.0, col.num),
#                          P_raw=rep(0.0, col.num), P_adj=rep(0.0, col.num))
tax_wilcox <- data.frame(Taxon=colnames(tax_dat), Median_Rhizo=rep(0.0, col.num),
                         Median_Soil=rep(0.0, col.num),
                         P_raw=rep(0.0, col.num), P_adj=rep(0.0, col.num))
# tax_wilcox[,c('Median_HZn', 'Median_LZn')] <- t(apply(tax_dat, 2, function(x) 
                                                      # tapply(x, group, median)))
tax_wilcox[,c('Median_Rhizo', 'Median_Soil')] <- t(apply(tax_dat, 2, function(x) 
                                                      tapply(x, group, median)))
tax_wilcox$P_raw <- apply(tax_dat, 2, function(x) 
                          wilcox.test(x~group, paired=FALSE, correct=FALSE)$p.value)
tax_wilcox <- na.omit(tax_wilcox)
tax_wilcox$P_adj <- p.adjust(tax_wilcox$P_raw, method="BH")
write.csv(tax_wilcox, "species_wilcox_median_YL_Rhizo_vs_Soil.csv", row.names=FALSE)
tax_dat$Group <- group
tax_dat <- gather(tax_dat, key="Taxon", value="Abund", -Group)
tax_dat$Taxon <- factor(x=tax_dat$Taxon, levels=unique(taxon), ordered=TRUE) 
#tax_wilcox <- tax_wilcox[tax_wilcox$P_raw<0.05,'Taxon']
#tax_dat <- tax_dat[tax_dat$Taxon %in% tax_wilcox,]
tax_dat.p <- ggplot(tax_dat, aes(x=Group, y=Abund)) +
  geom_boxplot(aes(fill=Group), color="black", width=0.6, size=0.3,
               outlier.size=0.3) +
  geom_signif(comparisons=list(c('HZn', 'LZn')), 
              test="wilcox.test", map_signif_level=TRUE, tip_length=0.02, 
              step_increase=0.05, textsize=2, size=0.3) +
  facet_wrap(~ Taxon, nrow=1, strip.position="top", scales="fixed") +
  labs(x=NULL, y="Relative abundance (%)", fill=NULL) +
  #scale_fill_manual(values=brewer.pal(12, "Paired")[c(4,8)]) +
  scale_fill_manual(values=brewer.pal(12, "Paired")[c(6,2)]) +
  #scale_y_continuous(limits=c(0,0.00015)) +
  theme(line=element_line(color="black", size=0.3), 
        text=element_text(color="black", size=8),
        axis.ticks=element_line(color="black", size=0.3),
        axis.text=element_text(color="black", size=8),
        panel.background=element_rect(fill="white", color="black"),
        panel.grid=element_blank(), legend.key=element_rect(fill="white"),
        legend.position="none")
ggsave("rhizo_OG_ZSB_selected_HLZn_wilcox_boxplot_YL.pdf", tax_dat.p, device="pdf",
       width=12, height=5, units="cm")
