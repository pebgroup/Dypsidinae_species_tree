wd <- "/Users/au265104/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/gis07.st.client.au.dk – SFTP/data_vol/wolf/Dypsis"
setwd(wd)

library(ape)
library(phytools)
library(classInt)
library(colorspace)
library(grDevices)
library(stringr)
library(adephylo)

source("/Users/au265104/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/gis07.st.client.au.dk – SFTP/home/au265104/scripts/dypsidinae/functions.R")
source("/Users/au265104/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/gis07.st.client.au.dk – SFTP/home/au265104/scripts/dypsidinae/load_trees.R")

# create astral tree with node labels corresponding to node numbers
astral_tree <- root(astral_tree, outgroup=c("1011", "1012"))
astral_tree <- ladderize(astral_tree, right=F)
astral_tree_nodeno <- astral_tree
#astral_tree_nodeno = root(astral_tree_nodeno, outgroup=c("1011", "1012"))
astral_tree_nodeno$node.label <- 1:length(astral_tree$node.label) + Ntip(astral_tree)
astral_tree_nodeno$edge.length[is.na(astral_tree_nodeno$edge.length)] <- 0.001
write.tree(astral_tree_nodeno, "final_tree_nofilter/astral/astral_tree_nodeno.tre")

#plot(astral_tree_nodeno, show.node.label=TRUE)

# create list of representation across gene trees: tipname, locusname

locusname <- c()
for(i in 1:length(gtnames)){
  locusname <- c(locusname, strsplit(gtnames[i],"_")[[1]][2])
}
rm(i)

data = matrix(nrow=0,ncol=2)
colnames(data) <- c("tipname","locusname")

for(i in 1:length(gts)){
  #i = 1
  for(tip in gts[[i]]$tip.label){
    data <- rbind(data,c(tip,locusname[i]))
    #print(tip)
  }
}

data <- as.data.frame(data)

write.table(data, file="final_tree_nofilter/astral/locusdata.csv", sep=",", col.names=T, row.names=F, quote=F)

# Run 
# export PYTHONPATH="${PYTHONPATH}:/home/au265104/scripts/physcripts"
# /data_vol/wolf/Dypsis/final_tree_nofilter/astral$ python ~/scripts/dypsidinae/calc_branchwise_decisiveness.py astral_tree_nodeno.tre locusdata.csv node_scores_outfile locus_scores_outfile branch_counts_outfile

node_scores <- read.table("final_tree_nofilter/astral/node_scores_outfile",sep=",")

classIntervals(node_scores$node_label, 6, "jenks") -> nodeclass
pal <- grDevices::rainbow(16)[11:16]
as.vector(findColours(nodeclass, pal)) -> node_scores$colour

astral_tree_node_scores <- astral_tree_nodeno
astral_tree_node_scores$node.label <- node_scores[as.character(astral_tree_nodeno$node.label),"node_label"]

plot(astral_tree_node_scores, show.node.label=TRUE)

edgecols <- rep("black",nrow(astral_tree_node_scores$edge))

# loop through internal nodes and colour branches
for(nl in as.character(astral_tree_nodeno$node.label)){
  if(nl %in% rownames(node_scores)){
    print(nl)
    edgecols[astral_tree_node_scores$edge[,2]==nl] <- node_scores[nl,"colour"]
  }
}

# loop through terminal nodes and colour branches
i = 1
for(nl in astral_tree_nodeno$tip.label){
  if(nl %in% rownames(node_scores)){
    print(nl)
    print(i)
    edgecols[astral_tree_node_scores$edge[,2]==i] <- node_scores[nl,"colour"]
    i = i +1
  }
}

# create tip colours to match branches
tipcols <- c()
i = 1
for(nl in astral_tree_nodeno$tip.label){
  if(nl %in% rownames(node_scores)){
    tipcols <- c(tipcols, node_scores[nl,"colour"])
    i = i +1
  }
}

#edgelty = rep(1,nrow(astral_tree_node_scores$edge))
#edgelty[astral_tree_node_scores$edge[,2]<=Ntip(astral_tree_node_scores)] <- 3

plot(astral_tree_for_figure, show.node.label=T)

astral_tree_for_figure$node.label -> supportvals
names(supportvals) <- 1:length(astral_tree_for_figure$node.label) + Ntip(astral_tree_for_figure)
edgelabs <- c()
for(i in 1:nrow(astral_tree_for_figure$edge)){
  dscndnt <- as.character(astral_tree_for_figure$edge[i,2])
  if(dscndnt %in% names(supportvals)){
    edgelabs <- c(edgelabs, supportvals[dscndnt])
  } else {
    edgelabs <- c(edgelabs, "")
  }
}

edgelabs <- rapply(as.list(as.numeric(edgelabs)), sprintf, fmt = "%0.2f", how = "replace")
edgelabs <- unlist(edgelabs)
edgelabs[edgelabs == "NA"] <- ""
edgelabs <- str_replace(edgelabs, "0.", ".")

astral_tree_for_figure <- astral_tree_node_scores
astral_tree_for_figure$tip.label = figurename_idx[astral_tree_for_figure$tip.label]
astral_tree_for_figure$node.label = astral_tree$node.label
astral_tree_for_figure$node.label[astral_tree_for_figure$node.label == "1"] <- ""
svg("/Users/au265104/Desktop/main_tree2.svg", height=15.3, width=8.25)
pdf("/Users/au265104/Desktop/main_tree2.pdf", height=15.3, width=8.25)
pdf("/Users/au265104/Desktop/main_tree2_broad.pdf", height=15.3, width=11)
#plot(astral_tree_for_figure, show.node.label=TRUE, cex = 0.4, align.tip.label = T, edge.color = edgecols)
plotwe(astral_tree_for_figure, cex = 0.45, align.tip.label = T, edge.color = edgecols, link.color = tipcols, label.offset = 0.08, x.lim = c(0,max(distRoot(astral_tree_for_figure))+1.3))
#nodelabels(text = astral_tree_for_figure$node.label, frame="none", cex=0.45, adj=c(-0.45,0.35))
edgelabels(text = edgelabs, frame="none", cex=0.45, adj=c(0.5,-0.30))
legend("topleft", legend=c("1-14","15-39","40-76","77-119","120-148","149-161"),fill=pal,title="No. gene trees")
add.scale.bar(x=0, y=20)
dev.off()




max(distRoot(astral_tree_for_figure))
