library(geiger)
library(ape)
library(phytools)

tree_optrimal <- read.tree("astral_tree_prelim_renamed.tre")
tree_optrimal$edge.length[is.na(tree_optrimal$edge.length)] <- 1
tree_optrimal <- root(tree_optrimal, node=getMRCA(tree_optrimal, c("Loxococcus-rupicola-SBL234-S35", "Loxococcus-rupicola-SBL8-S7")))

tree_gt0.5 <- read.tree("astral_tree_prelim_gt0.5_renamed.tre")
tree_gt0.5$edge.length[is.na(tree_gt0.5$edge.length)] <- 1
tree_gt0.5 <- root(tree_gt0.5, node=getMRCA(tree_gt0.5, c("Loxococcus-rupicola-SBL234-S35", "Loxococcus-rupicola-SBL8-S7")))

association <- cbind(tree_optrimal$tip.label, tree_optrimal$tip.label)

coplot <- cophylo(tree_optrimal, tree_gt0.5, assoc=association, rotate=TRUE)

png(filename="cophyloplot2.png", width=3500, height=3500)

cophyloplot(coplot[[1]][[1]], coplot[[1]][[2]], assoc=coplot[[2]], use.edge.length = FALSE, space = 200,
            length.line = 30, gap = 3, type = "phylogram", rotate = FALSE,
            col = par("fg"), lwd = par("lwd"), lty = par("lty"),
            show.tip.label = TRUE, cex=.5)  #rotate lets you rotate in plot panel

dev.off()

#Check tree characters; is it rooted what are the tip names etc.               
SaxifragaTree   
str(SaxifragaTree)
SaxifragaTree$tip.label # shows vector tip.label; tip names and tip numbers
SaxifragaTree$Nnode
plot(SaxifragaTree)

#do we need to properly re-root the tree?
findMRCA(SaxifragaTree,c("Itea_coriacea_S372", "Itea_chinensis_S371")) #find ancestral node: 615
SaxifragaTreeRoot <- reroot(SaxifragaTree, 615, position=NULL, interactive=FALSE)

findMRCA(SaxifragaTreeTRIM,c("Itea_coriacea_S372", "Itea_chinensis_S371"))#find ancestral node: 508
SaxifragaTreeTRIMRoot <- reroot(SaxifragaTreeTRIM, 508, position=NULL, interactive=FALSE)
plot(SaxifragaTreeTRIMRoot)



#creation of the association matrix:
#what does it need to look like? 
#Any tutorial makes us set the same trees against each other, not the different trees
#Either way the association matrix does not our absence of cross-links
association <- cbind(SaxifragaTreeRoot$tip.label, SaxifragaTreeRoot$tip.label)


#make cophyloplot
png(filename="cophyloplot.png")
cophyloplot(SaxifragaTreeRoot, SaxifragaTreeTRIMRoot, assoc = association, use.edge.length = FALSE, space = 0,
            length.line = 30, gap = 3, type = "phylogram", rotate = FALSE,
            col = par("fg"), lwd = par("lwd"), lty = par("lty"),
            show.tip.label = FALSE)  #rotate lets you rotate in plot panel
dev.off()


####  Result: the trees have the same tips, same root, but the combination of association matrix and cophyloplot do not result in links between the same tip labels!!




##second option: cophylo function from phytoolS, but it does not plot a figure. Can be run without cophyloplot to no avail (finds copies of names automatically)
coplot <- cophylo(SaxifragaTreeRoot, SaxifragaTreeTRIMRoot, assoc=association, rotate=TRUE)
png(filename="cophyloplot2.png")
plot(coplot)
dev.off()
