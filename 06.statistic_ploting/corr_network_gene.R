#clean
rm(list=ls())
gc()
setwd(dir=getwd())

#load package
library(psych)
library(RColorBrewer)
library(igraph)
library(ggplot2)

#read data
tax <- read.csv("combined_42.species_abundance", header=T, sep="\t", row.names=1)
tax <- as.data.frame(t(tax))
HZn_vs_LZn <- read.csv("species_wilcox_YS_HZn_vs_LZn.csv", header=T)
HZn_vs_Soil <- read.csv("species_wilcox_YS_HZn_vs_Soil.csv", header=T)
LZn_vs_Soil <- read.csv("species_wilcox_YS_LZn_vs_Soil.csv", header=T)

#correlation among taxa
Rhizo <- intersect(HZn_vs_Soil[HZn_vs_Soil$Diff=="HZn",'Taxon'],
                   LZn_vs_Soil[LZn_vs_Soil$Diff=="LZn",'Taxon'])
HZn_Enr <- intersect(Rhizo, HZn_vs_LZn[HZn_vs_LZn$Diff=="HZn",'Taxon'])
LZn_Enr <- intersect(Rhizo, HZn_vs_LZn[HZn_vs_LZn$Diff=="LZn",'Taxon'])
tax.dat <- tax[grepl("YS.[1-4].[^B]", rownames(tax)), colnames(tax) %in% Rhizo]
#cor.mx <- corr.test(x=tax.dat, method="spearman", adjust="BH", ci=FALSE)
#cor.mx$r[cor.mx$p>0.05] <- 0
#write.csv(cor.mx$r, "species_corr_matrix_HLZn_YL.csv", row.names=T)
cor.mx <- cor(tax.dat, method="spearman")
cor.mx[abs(cor.mx)<0.7] <- 0
write.csv(cor.mx, "species_corr_matrix_HLZn_YS.csv", row.names=T)
cor.dat <- data.frame(Taxon=Rhizo, Group=rep("", ncol(tax.dat)),
                      Median=rep(0.0, ncol(tax.dat)))
cor.dat[cor.dat$Taxon %in% HZn_Enr, 'Group'] <- "HZn"
cor.dat[cor.dat$Taxon %in% LZn_Enr, 'Group'] <- "LZn"
cor.dat$Median <- apply(tax.dat, 2, median)
write.csv(cor.dat, "species_corr_abund_HLZn_YS.csv", row.names=F)

#draw co-occurrence network
cor.mx <- read.csv("genus_corr_matrix_HLZn_YL.csv", header=T, row.names=1)
cor.dat <- read.csv("genus_corr_abund_HLZn_YL.csv", header=T, row.names=1)
colnames(cor.mx) <- rownames(cor.mx)
cor.mx <- as.matrix(cor.mx)
cor.mx[abs(cor.mx)<0.6] <- 0
net <- graph_from_adjacency_matrix(cor.mx, mode="undirected", weighted=T, diag=F)
small <- rownames(subset(cor.dat, Median<0.01))
net <- delete.vertices(net, V(net)$name %in% small) #omit low-abund taxa
net <- delete.vertices(net, V(net)[degree(net)==0]) #omit taxa with no link
net.wt <- E(net)$weight #weight as Spearman correlation
E(net)$color = as.character(ifelse(net.wt>0, brewer.pal(12, "Paired")[9], 
                                   brewer.pal(12, "Paired")[1]))
E(net)$width = abs(net.wt)
net.sz <- cor.dat[V(net)$name, 'Median']
V(net)$size = log2(net.sz*64+2)
net.gp <- as.factor(cor.dat[V(net)$name, 'Group'])
levels(net.gp) = c(brewer.pal(8, "Set2")[c(8,2,1)])
V(net)$color = as.character(net.gp)
#color vertex by module
#fc = cluster_fast_greedy(net, weights=NULL)
#modularity = modularity(net, membership(fc))
#comps = membership(fc)
#colbar = brewer.pal(max(comps), "Set2")
#V(net)$color = colbar[comps] 
pdf("species_network_HLZn_YL.pdf", width=12, height=12, family="Helvetica", pointsize=10)
E(net)$weight=NA
set.seed(2020)
plot(net, vertex.frame.color="black", vertex.label=NA, margin=c(0,0,0,0),
     edge.curved=FALSE, layout=layout.fruchterman.reingold)
dev.off()

#statistics of co-occurrence network
#connectance indicates the ratio of actual to total potential edges
edge_density(net, loops=F)
#degree per vertex
head(sort(degree(net), decreasing=TRUE)) #top 10
mean(degree(net)) #average
#average path length between two vertexes
average.path.length(net)
#diameter of the network
diameter(net, directed=F, unconnected=T, weights=NULL)
#edge connectivity or group adhesion
edge_connectivity(net)
#clustering coefficient or transitivity, the network is clustered or dispersed
transitivity(net)
no.clusters(net)
#betweeness and degree centralizality
head(sort(centralization.betweenness(net)$res, decreasing=TRUE)) #top 10
centralization.betweenness(net)$centralization #average
head(sort(centralization.degree(net)$res, decreasing=TRUE)) #top 10
centralization.degree(net)$centralization #average

#statistics of degrees and edges by group
length(E(net)) #number of edges
length(V(net)) #number of vertexes
sum(net.wt>0) #number of positive correlations
sum(net.wt<0) #number of negative correlations
head(sort(degree(net), decreasing=TRUE)) #top 10
mean(degree(net)) #average
dat <- as_data_frame(net)
HZn <- rownames(subset(cor.dat, Group=="HZn"))
LZn <- rownames(subset(cor.dat, Group=="LZn"))
Rhizo <- rownames(subset(cor.dat, Group==""))
dat.eg <- data.frame(Link=rep(c("Positive","Negative"), 6), 
                     Group=rep(c("HZn","LZn","H_LZn","Rhizo","HZn_Rhizo","LZn_Rhizo"),
                               each=2),
                     Number=rep(0, 12))
dat.eg$Number[1] <- nrow(subset(dat, from %in% HZn & to %in% HZn & weight > 0))
dat.eg$Number[2] <- nrow(subset(dat, from %in% HZn & to %in% HZn & weight < 0))
dat.eg$Number[3] <- nrow(subset(dat, from %in% LZn & to %in% LZn & weight > 0))
dat.eg$Number[4] <- nrow(subset(dat, from %in% LZn & to %in% LZn & weight < 0))
dat.eg$Number[5] <- nrow(subset(dat, (from %in% HZn & to %in% LZn |
                                        from %in% LZn & to %in% HZn) & weight > 0))
dat.eg$Number[6] <- nrow(subset(dat, (from %in% HZn & to %in% LZn |
                                        from %in% LZn & to %in% HZn) & weight < 0))
dat.eg$Number[7] <- nrow(subset(dat, from %in% Rhizo & to %in% Rhizo & weight > 0))
dat.eg$Number[8] <- nrow(subset(dat, from %in% Rhizo & to %in% Rhizo & weight < 0))
dat.eg$Number[9] <- nrow(subset(dat, (from %in% HZn & to %in% Rhizo |
                                        from %in% Rhizo & to %in% HZn)& weight > 0))
dat.eg$Number[10] <- nrow(subset(dat, (from %in% HZn & to %in% Rhizo |
                                         from %in% Rhizo & to %in% HZn)& weight < 0))
dat.eg$Number[11] <- nrow(subset(dat, (from %in% LZn & to %in% Rhizo |
                                         from %in% Rhizo & to %in% LZn)& weight > 0))
dat.eg$Number[12] <- nrow(subset(dat, (from %in% LZn & to %in% Rhizo |
                                         from %in% Rhizo & to %in% LZn)& weight < 0))
dat.eg.p <- ggplot(dat.eg, aes(x=Group, y=Number)) +
  geom_bar(aes(fill=Link), stat="identity", position="stack", 
           color="black", width=0.6, size=0.3) +
  scale_fill_manual(values=brewer.pal(12, "Paired")[c(1,9)]) +
  scale_x_discrete(limits=c("HZn","LZn","H_LZn","Rhizo","HZn_Rhizo","LZn_Rhizo")) +
  labs(x=NULL, y="Number of edges", fill="Correlation") +
  theme(line=element_line(color="black", size=0.3), 
        text=element_text(color="black", size=8),
        axis.ticks=element_line(color="black", size=0.3),
        axis.text=element_text(color="black", size=8),
        panel.background=element_rect(fill="white", color="black"),
        panel.grid=element_blank(), legend.key=element_rect(fill="white"))
ggsave("genus_network_edge_barplot_YL.pdf", dat.eg.p, width=14, height=8, units="cm")
