#clean
rm(list=ls())
gc()
setwd("D:/项目试验方案数据/小麦微生物组/指标数据/宏基因组分析结果/06.function")

#load package
library(venn)
library(ggplot2)
library(ggtern)
library(RColorBrewer)

#venn diagram of enriched taxa
HZn_vs_LZn <- read.csv("KO_wilcox_YL_HZn_vs_LZn.csv", header=T)
HZn_vs_Soil <- read.csv("KO_wilcox_YL_HZn_vs_Soil.csv", header=T)
LZn_vs_Soil <- read.csv("KO_wilcox_YL_LZn_vs_Soil.csv", header=T)
pdf(file="species_venn_rhizo_vs_HZn_vs_LZn_YL.pdf", width=7, height=7, pointsize=16)
venn(list('HZn>Soil'=HZn_vs_Soil[HZn_vs_Soil$Group=="HZn", 'Taxon'],
          'HZn>LZn'=HZn_vs_LZn[HZn_vs_LZn$Group=="HZn", 'Taxon'], 
          'LZn>Soil'=LZn_vs_Soil[LZn_vs_Soil$Group=="LZn", 'Taxon']),
     zcolor=brewer.pal(12, "Paired")[c(6,2,4)], plotsize=15, ilcs=1, sncs=1.5, box=FALSE)
dev.off()

YL_rhizo <- read.csv("OG_wilcox_YL_Rhizo_vs_Soil.csv", header=T)
#YL_rhizo <- read.csv("rhizo_OG_wilcox_YL_HZn_vs_LZn.csv", header=T)
#YL_HLZn <- read.csv("KO_wilcox_YL_HZn_vs_LZn.csv", header=T)
#YS_rhizo <- read.csv("KO_wilcox_YS_Rhizo_vs_Soil.csv", header=T)
#YS_HLZn <- read.csv("KO_wilcox_YS_HZn_vs_LZn.csv", header=T)
#pdf(file="KO_venn_rhizo_vs_soil_YL.pdf", width=7, height=7, pointsize=16)
#venn(list('Rhizo'=YL_rhizo[YL_rhizo$Group=="Rhizo", 'Taxon'],
#          'Zother'=YL_rhizo[YL_rhizo$Group=="",'Taxon'],
#          'Soil'=YL_rhizo[YL_rhizo$Group=="Soil", 'Taxon']),
#     zcolor=brewer.pal(12, "Paired")[c(6,2,4)], plotsize=15, ilcs=1, sncs=1.5, box=FALSE)
#dev.off()
number <- c(length(YL_rhizo[YL_rhizo$Group=='Rhizo', 'Taxon']),
            length(YL_rhizo[YL_rhizo$Group=='', 'Taxon']),
            length(YL_rhizo[YL_rhizo$Group=='Soil', 'Taxon']))
#number <- c(length(YL_rhizo[YL_rhizo$Group=='HZn', 'Taxon']),
#            length(YL_rhizo[YL_rhizo$Group=='', 'Taxon']),
#            length(YL_rhizo[YL_rhizo$Group=='LZn', 'Taxon']))
perc <- paste(round(number/sum(number) * 100, 2), "%")
group <- c('Rhizo', 'Not_enrich', 'Soil')
#group <- c('HZn', 'Not_enrich', 'LZn')
cols <- c(brewer.pal(12, "Paired")[4], 'grey90', brewer.pal(12, "Paired")[8])
#cols <- c(brewer.pal(12, "Paired")[6], 'grey90', brewer.pal(12, "Paired")[4])
pdf(file="OG_pie_Rhizo_vs_Soil_YL.pdf", width=7, height=7, pointsize=16)
#pdf(file="rhizo_OG_pie_HZn_vs_LZn_YL.pdf", width=7, height=7, pointsize=16)
pie(number, labels=paste(group, number, perc, sep=" "), col=cols, 
    main="Number of enriched OGs", clockwise=TRUE)
dev.off()

#ternary plot of enriched taxa
tax <- read.csv("combined_42.KEGG_KO_abundance", header=T, sep="\t", row.names=1)
tax <- tax[,grepl("YL", colnames(tax))]
group <- as.factor(rep(c('LZn', 'HZn', 'LZn', 'LZn', 'HZn', 'HZn', 'Soil'), 4))
tax_dat <- data.frame(Taxon=rownames(tax), HZn=rep(0.0, nrow(tax)), 
                      LZn=rep(0.0, nrow(tax)), Soil=rep(0.0, nrow(tax)), 
                      Median=rep(0.0, nrow(tax)), Group=rep('', nrow(tax)))
tax_dat[,2:4] <- t(apply(tax, 1, function(x) tapply(x, group, median)))
tax_dat[,5] <- apply(tax, 1, median)
Rhizo <- intersect(HZn_vs_Soil[HZn_vs_Soil$Group=="HZn",'Taxon'],
                   LZn_vs_Soil[LZn_vs_Soil$Group=="LZn",'Taxon'])
HZn_Enr <- intersect(Rhizo, HZn_vs_LZn[HZn_vs_LZn$Group=="HZn",'Taxon'])
LZn_Enr <- intersect(Rhizo, HZn_vs_LZn[HZn_vs_LZn$Group=="LZn",'Taxon'])
Soil_Enr <- intersect(HZn_vs_Soil[HZn_vs_Soil$Group=="Soil",'Taxon'],
                   LZn_vs_Soil[LZn_vs_Soil$Group=="Soil",'Taxon'])
tax_dat[tax_dat$Taxon %in% HZn_Enr, 'Group'] <- 'HZn'
tax_dat[tax_dat$Taxon %in% LZn_Enr, 'Group'] <- 'LZn'
tax_dat[tax_dat$Taxon %in% Soil_Enr, 'Group'] <- 'Soil'
tax.p <- ggtern(tax_dat, aes(x=HZn, y=Soil, z=LZn)) +
  geom_point(aes(color=Group, size=Median)) +
  scale_color_manual(values=c("grey90", brewer.pal(12, "Paired")[c(6,2,4)])) +
  labs(color="enriched", size="abundance") +
  theme_bw() +
  theme(line=element_line(color="black", size=0.3), 
        text=element_text(color="black", size=8),
        axis.line=element_line(color="black", size=0.3),
        axis.ticks=element_line(color="black", size=0.3),
        axis.text=element_text(color="black", size=8),
        legend.key=element_rect(fill="white"))
ggsave("KO_thernary_enriched_YL.pdf", tax.p, device="pdf", 
       width=18, height=12, units="cm")
