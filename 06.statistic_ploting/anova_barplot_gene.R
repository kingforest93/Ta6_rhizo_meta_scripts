#clean
rm(list=ls())
gc()
setwd(dir=getwd())

#load package
library(tidyr)
library(ggplot2)
library(agricolae)
library(RColorBrewer)

#read tax count matrix
tax <- read.csv("combined_42.eggNOG_OG_abundance", sep="\t", header=TRUE, row.names=1)
list <- read.table("clipboard", header=FALSE)

#anova of each taxon by cultivar
Diff <- list[,1]
tax_res <- tax[rownames(tax) %in% Diff, grepl("YS", colnames(tax))]
total <- apply(tax_res, 2, sum)
tax_res <- rbind(tax_res, total)
rownames(tax_res)[nrow(tax_res)] <- 'Total'
tax_res <- tax_res[order(apply(tax_res, 1, function(x)max(x, na.rm=TRUE)), 
                         decreasing=TRUE),]
tax_res <- as.data.frame(t(tax_res))
anno_df <- data.frame(Taxon=rep(colnames(tax_res), each=7),
                      Group=rep(c('LZn','HZn','LZn','LZn','HZn','HZn','Soil'), 
                                ncol(tax_res)),
                      Cult=rep(rep("", 7), ncol(tax_res)),
                      Mean=rep(rep(0, 7), ncol(tax_res)),
                      SE=rep(rep(0, 7), ncol(tax_res)),
                      Sig=rep(rep("", 7), ncol(tax_res)))
tax_res$Cultivar <- as.factor(rep(c('10_yannong0428','14_zhoumai24','17_hengguan35', 
                                    '18_jinan17','3_bei9','4_xinong3517','Soil'), 4))
tax_res$Block <- as.factor(rep(1:4, each=7))
for (i in 1:(ncol(tax_res)-2)) {
  tem <- aov(tax_res[,i] ~ Cultivar + Block, data=tax_res)
  cmp <- duncan.test(tem, 'Cultivar', alpha=0.05, group=TRUE)
  anno_df[(7*(i-1)+1):(7*i), 'Cult'] <- rownames(cmp$means)
  anno_df[(7*(i-1)+1):(7*i), 'Mean'] <- cmp$means[,1]
  anno_df[(7*(i-1)+1):(7*i), 'SE'] <- cmp$means$std / sqrt(cmp$means$r)
  anno_df[(7*(i-1)+1):(7*i), 'Sig'] <- cmp$groups[rownames(cmp$means), 'groups']
}
anno_df$Group <- as.factor(anno_df$Group)
anno_df$Cult <- as.factor(anno_df$Cult)
anno_df.p <- ggplot(anno_df, aes(x=Cult, y=Mean, ymin=Mean-SE, ymax=Mean+SE)) +
  geom_col(aes(fill=Group), color="black", width=0.6, size=0.3) +
  geom_errorbar(color="black", width=0.15, size=0.3) +
  geom_text(aes(y=Mean+SE*2, label=Sig), size=2) +
  facet_wrap(~ Taxon, nrow=21, strip.position="top", scales="free") +
  labs(x="Cultivar", y="Relative abundance (%)", fill=NULL) +
  scale_fill_manual(values=brewer.pal(12, "Paired")[c(6,2,4)]) +
  theme(line=element_line(color="black", size=0.3), 
        text=element_text(color="black", size=8),
        axis.ticks=element_line(color="black", size=0.3),
        axis.text=element_text(color="black", size=8),
        axis.text.x=element_text(color="black", size=8, angle=45, hjust=1),
        panel.background=element_rect(fill="white", color="black"),
        panel.grid=element_blank(), legend.key=element_rect(fill="white"),
        legend.position="top")
ggsave("siderophore_eggNOG_OG_anova_boxplot_cultivar_YS.pdf", anno_df.p, device="pdf",
       width=100, height=125, units="cm")
