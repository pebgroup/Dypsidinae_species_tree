setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(ape)
library(phytools)
library(gtools)
library(ggtree)
library(cowplot)
library(ggimage)
library(treeio)
# you probably don't need all the above packages but you certainly need ggimage, which requires an up-to-date R

tree <- read.tree("trees/astral_tree_prelim_t2_renamed.tre")
tree <- root(tree, node=getMRCA(tree, c("Loxococcus-rupicola-SBL234-S35", "Loxococcus-rupicola-SBL8-S7")), edgelabel = TRUE)
write.tree(tree, "trees/astral_tree_prelim_t2_renamed_rerooted.tre")
ape::write.nexus(tree, file = "test.nex", translate=FALSE)

########### Example plot pie charts with Q scores on ASTRAL topologies



# import annotated ASTRAL tree (I like to use a tree rooted with pxrr directly after getting it from ASTRAL, before any editing)
# Then the tree has to be formated in nexus (for instance open it in Figtree and export it as nexus including annotations - don't do any other change in figtree, just open and export)
# Then it has to be edited this way (for example using Notepad++):
# replace [&label="[ by [&
# replace ]"] by ]
# replace =([\d\.]+)\- by =\1
# replace =([\d\.]+E\-\d+)\- by =\1
# Then you can import the edited nexus tree using the read.beast function:

AS_all <- read.beast("trees/astral_tree_prelim_t2_renamed_rerooted_ed.nex")

# plot and ladderize the tree, without using the ASTRAL branch lengths
p <- ggtree(AS_all@phylo, ladderize=T, branch.length = "none") + geom_tiplab(size=2.5, hjust= -0.05) + xlim_tree(70) + ggtitle("ASTRAL optrimal") 

# check node labels if necessary
#p + geom_text2(aes(subset=!isTip, label=node), hjust=-.3, size = 4)

# format data to make pie charts (or any plot) of the Q scores (could do it for any other data from AS_all@data)
Q1 <- as.numeric(AS_all@data$q1) * 100
Q <- as.data.frame(Q1)
Q$Q2 <- as.numeric(AS_all@data$q2) * 100
Q$Q3 <- as.numeric(AS_all@data$q3) * 100
Q$node <- AS_all@data$node

# make barplots
# bars <- nodebar(Q, cols=1:3, position='dodge', color=c(Q1='red', Q2='cyan', Q3='gray'))
# inset(p, bars, x='node', width=8, height=2)

# make pie charts
pies <- nodepie(Q, cols=1:3, color=c(Q1='blue', Q2='green', Q3='red')) #, alpha=.6) # change style as needed, alpha is for transparency

png(filename="plots/piecharts.png", width=1750, height=3000)

inset(p, pies, width=0.1, height=0.1)#, hjust=-.6)

dev.off()

# pies give the quartet support: percentage of quartets agreeing: 
# with the branch (blue), 
# with the second alternative RS|LO (green), 
# and with the last alternative RO|LS (red).




# FROM HERE SIDONIE

##### plot QS to see their values
AS_A <- read.beast("A/A_SpeciesTree_BP10_annotQ_rooted2_ed.nex")

Q1 <- as.numeric(AS_A@data$q1) * 100
Q <- as.data.frame(Q1)
Q$Q2 <- as.numeric(AS_A@data$q2) * 100
Q$Q3 <- as.numeric(AS_A@data$q3) * 100
Q$node <- AS_A@data$node
node <- AS_A@data$node
allQS <- as.data.frame(node)
for (x in (1:length(allQS$node))) {allQS$QS[x] <- paste(round(Q$Q1[x]), "_", round(Q$Q2[x]), "_", round(Q$Q3[x]), sep = "")}

pP <- ggtree(AS_A@phylo, ladderize=T, branch.length = "none") %<+% allQS + 
  geom_nodelab(aes(x=branch, label=QS), vjust=-0.5, size=3) +
  geom_tiplab(size=3, hjust= -0.05) + xlim_tree(70) + 
  ggtitle("ASTRAL, plastid regions")

pdf("ASTRAL_BP10_QSvalues.pdf", 15, 15)
plot_grid(pP, pN, ncol=2)
dev.off()