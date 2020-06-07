setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(geiger)
library(phytools)
library(colorspace)
library(ggplot2)
#library(BiocManager)
library(vctrs)
library(ggtree)

# Load trees
tree_optrimal <- read.tree("trees/astral_tree_prelim_renamed.tre")
tree_optrimal$edge.length[is.na(tree_optrimal$edge.length)] <- 1
tree_optrimal <- ape::root(tree_optrimal, node=getMRCA(tree_optrimal, c("Loxococcus-rupicola-SBL234-S35", "Loxococcus-rupicola-SBL8-S7")), edgelabel = TRUE)

tree_gt0.5 <- read.tree("trees/astral_tree_prelim_gt0.5_renamed.tre")
tree_gt0.5$edge.length[is.na(tree_gt0.5$edge.length)] <- 1
tree_gt0.5 <- ape::root(tree_gt0.5, node=getMRCA(tree_gt0.5, c("Loxococcus-rupicola-SBL234-S35", "Loxococcus-rupicola-SBL8-S7")), edgelabel = TRUE)

# Load gene tree stats (number of gene trees per species) and PoM groups
gtstats <- read.csv("genetrees_per_taxon_prelim.csv",sep=";")
gtstats <- gtstats[!is.na(gtstats$No.of.genes) & gtstats$Tip.label %in% tree_optrimal$tip.label,]
nogt <- gtstats$No.of.genes
names(nogt) <- gtstats$Tip.label

# Define twelve categories for colour scale
breaks <- seq(from=min(nogt), to=max(nogt), length.out = 13)
breaks[1] <- -Inf
breaks[13]  <- Inf

nogt_class <- cut(nogt, breaks = breaks, labels=1:12)

# Define colours
nogt_col <- diverge_hcl(12,h=c(255,330),l=c(40,90))[as.numeric(as.vector(nogt_class))]
names(nogt_col) <- gtstats$Tip.label


# Plot cophyloplot
association <- cbind(tree_optrimal$tip.label, tree_optrimal$tip.label)

coplot <- cophylo(tree_optrimal, tree_gt0.5, assoc=association, rotate=TRUE)

png(filename="plots/cophyloplot3.png", width=3500, height=3500)

cophyloplot(coplot[[1]][[1]], coplot[[1]][[2]], assoc=coplot[[2]], use.edge.length = FALSE, space = 200,
            length.line = 30, gap = 3, type = "phylogram", rotate = FALSE,
            col=nogt_col[association[,1]], lwd = par("lwd"), lty = par("lty"),
            show.tip.label = TRUE, cex=.5)  #rotate lets you rotate in plot panel

breaks_for_legend_floor <- floor(seq(from=min(nogt), to=max(nogt), length.out = 13))
breaks_for_legend_ceiling <- ceiling(seq(from=min(nogt), to=max(nogt), length.out = 13))
leg <- c()
for(i in 1:12){
  leg <- c(leg, paste(breaks_for_legend_ceiling[i],"-",breaks_for_legend_floor[i+1])) 
}

legend("top", legend=leg, fill=diverge_hcl(12,h=c(255,330),l=c(40,90)), cex=2, horiz=TRUE)

dev.off()

# Phylo barplot

png(filename="plots/phylobarplot.png", width=1500, height=2000)

plotTree.barplot(tree_optrimal, nogt, scale=1000, width=10)

dev.off()



# PoM groups 

# Create PoM group presence/absence matrix
pomgroup <- gtstats$PoM.Group
names(pomgroup) <- gtstats$Tip.label
pomgroup <- pomgroup[!is.na(pomgroup)]
heatmapData <- matrix(nrow = length(tree_optrimal$tip.label), ncol=18)
rownames(heatmapData) <- tree_optrimal$tip.label
heatmapData[,] <- 0
for(i in 1:length(pomgroup)){
  heatmapData[names(pomgroup[i]), pomgroup[i]] <- 1
}
rn <- rownames(heatmapData)
heatmapData <- as.data.frame(heatmapData)
heatmapData <- as.data.frame(sapply(heatmapData, as.character))
rownames(heatmapData) <- rn
colnames(heatmapData) <- paste("G",1:18,sep="")

png(filename="plots/pomgroups2.png", width=3000, height=1750)

p <- ggtree(tree_optrimal, layout="rectangular") + geom_tiplab(size=3, align=TRUE, linetype='dashed', linesize=.3) + geom_text2(aes(subset = !isTip & as.numeric(label) < 1, label=label, color="red"))
gheatmap(p, heatmapData, offset=1.5, width=.4, font.size=3) +
  scale_fill_manual(breaks=c("0", "1"), values=c("lightgrey", "black"), name="PoM group")

dev.off()