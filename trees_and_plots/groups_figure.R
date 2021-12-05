library(phytools)
library(geiger)
library(stringr)

# Assuming script location within Git repo as working directory
source("functions.R")
source("load_trees.R")

traitData3 <- read.csv2("~/OneDrive - Aarhus Universitet/PROJECTS/65 Dypsis systematics paper/~manuscript/figures/trait figure/DypsidinaeTraitData.csv",header = TRUE)
traitData3 <- traitData3[!traitData3$SpecName=="",]
traitData3$SpecName[traitData3$SpecName == "Dypsis Leucomalla"] <- "Dypsis leucomalla"
traitData3 <- traitData3[!traitData3$SpecName=="Loxococcus rupicola",] #exclude Loxococcus

# manually ladderize astral_tree_to_ladderize.tre in FigTree and save as astral_tree_lad.tre (ladderizing in R messed up the plot)
# tree <- astral_tree
# tree$edge.length[is.na(tree$edge.length)] <- 0.001
# write.tree(tree,paste(data_dir, "/final_tree_nofilter/astral/astral_tree_to_ladderize.tre", sep="")) 

tree <- read.tree(paste(data_dir, "/final_tree_nofilter/astral/astral_tree_lad.tre", sep=""))
# remove outgroup and redundant samples (where more than one individual per species)
tree <- drop.tip(tree, c("0194", "0196", "0199", "0202", "0204", "1012", "1011"))

# load PoM group data
groups <- read.table(paste(data_dir, "/POM_groups.csv", sep=""), sep=";", colClasses = c("character", "character", "numeric"))
groups <- groups[1:189,]
groups <- groups[groups$V2!="",]
rownames(groups) <- groups$V2

# rename tree with figurenames
figurename <- figurename[1:174,]
rownames(figurename) <- figurename$V1
tree$tip.label <- figurename[tree$tip.label,"V2"]
# name.check(tree, groups)

# Build plot
############

groups <- groups[groups$V1 %in% figurename$V1,]
secapr2figname <- figurename$V2
names(secapr2figname) <- figurename$V1
#cbind(secapr2figname[groups$V1], rownames(groups))
rownames(groups) <- secapr2figname[groups$V1]

char <- groups$V3
names(char) <- rownames(groups)
char <- char[tree$tip.label]

#figurepath <- "/Users/au265104/OneDrive - Aarhus Universitet/PROJECTS/65 Dypsis systematics paper/~manuscript/figures/group figure"
figurepath <- "/Users/au265104/OneDrive - Aarhus Universitet/PROJECTS/65 Dypsis systematics paper/~manuscript/figures/group figure/new"
  
pdf(paste(figurepath, "/groupplot_new.pdf", sep=""), width=8.3, height = 11.7)

n <- length(tree$tip.label) # number of tips in the tree
plot.phylo(tree, label.offset=5.5, align.tip.label = T, cex = 0.45)

xoffset = 9.15
increment = .3
scf = 1.8

for(i in 1:18){
  #Group i
  grid <- rep(".", n)
  grid[is.na(char)] <- NA
  points(x=rep(xoffset, n), y=1:n, pch=grid, cex=.6, col = "grey")
  symbol <- rep(NA, n)
  symbol[char==i] <- 19
  #symbol[is.na(char)] <- NA
  points(x=rep(xoffset, n), y=1:n, pch=symbol, cex=.6)
  xoffset <- xoffset + increment
}

dev.off()

