wd <- "/Users/au265104/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/gis07.st.client.au.dk – SFTP/data_vol/wolf/Dypsis"
setwd(wd)

figurepath <- "/Users/au265104/OneDrive - Aarhus Universitet/PROJECTS/65 Dypsis systematics paper/~manuscript/figures/group figure"

library(ape)
library(phytools)
library(classInt)
library(colorspace)
library(grDevices)
library(stringr)
library(adephylo)
library(ggtree)
library(picante)

source("/Users/au265104/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/gis07.st.client.au.dk – SFTP/home/au265104/scripts/dypsidinae/functions.R")
source("/Users/au265104/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/gis07.st.client.au.dk – SFTP/home/au265104/scripts/dypsidinae/load_trees.R")

astral_tree <- root(astral_tree, outgroup=c("1011", "1012"))
astral_tree <- ladderize(astral_tree, right=F)

astral_tree_for_figure <- astral_tree
astral_tree_for_figure$tip.label = figurename_idx[astral_tree_for_figure$tip.label]

groups <- read.table("POM_groups.csv", sep=";", colClasses = c("character", "character", "numeric"))

# Create PoM group presence/absence matrix
pomgroup <- groups$V3
names(pomgroup) <- groups$V2
pomgroup <- pomgroup[!is.na(pomgroup)]
pomgroup <- pomgroup[names(pomgroup) %in% astral_tree_for_figure$tip.label]

heatmapData <- matrix(nrow = length(astral_tree_for_figure$tip.label), ncol=18)
rownames(heatmapData) <- astral_tree_for_figure$tip.label
heatmapData[,] <- 0
for(i in 1:length(pomgroup)){
  heatmapData[names(pomgroup[i]), pomgroup[i]] <- 1
}
rn <- rownames(heatmapData)
heatmapData <- as.data.frame(heatmapData)
heatmapData <- as.data.frame(sapply(heatmapData, as.character))
rownames(heatmapData) <- rn
colnames(heatmapData) <- paste("G",1:18,sep="")

svglite(filename="/Users/au265104/OneDrive - Aarhus Universitet/PROJECTS/65 Dypsis systematics paper/~manuscript/figures/group figure/pomgroup.svg", width=3000, height=1750)

p <- ggtree(astral_tree_for_figure, layout="rectangular") + geom_tiplab(size=3, align=TRUE, linetype='dashed', linesize=.3)# + geom_text2(aes(subset = !isTip & as.numeric(label) < 1, label=label, color="red"))
gheatmap(p, heatmapData, offset=1.5, width=.4, font.size=3) +
  scale_fill_manual(breaks=c("0", "1"), values=c("lightgrey", "black"), name="PoM group")

dev.off()

plot(astral_tree_for_figure)
extract.clade(astral_tree_for_figure, getMRCA(astral_tree_for_figure, c("D. laevis", "V. dransfieldii"))) -> astral_tree_for_ps
