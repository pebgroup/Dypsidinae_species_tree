#!/usr/bin/Rscript

library(ape)

alnstats <- read.csv("alignments_trimmed/alnstats.csv", sep=";")

Nnode <- rep(NA,nrow(alnstats))
Nnode70 <- rep(NA,nrow(alnstats))

for(i in 1:nrow(alnstats)){
	name <- alnstats[i, "Alignment_name"]
	treepath <- paste("raxmltrees/", name, ".raxml.support", sep="")
	tree <- read.tree(treepath)
	Nnode[i] <- di2multi(tree, tol=0.0000011)$Nnode
	Nnode70[i] <- sum(as.numeric(tree$node.label)>=70, na.rm=T)
	}

alnstats$Nnode <- Nnode
alnstats$Nnode70 <- Nnode70

rm(Nnode, Nnode70)
attach(alnstats)

Nnode_perc <- Nnode/(No_of_taxa-2)
Nnode70_perc <- Nnode70/(No_of_taxa-2)

jpeg(filename="stats/Nnodes_hist.jpeg")
hist(Nnode_perc, main = "percentage of resolved nodes")
dev.off()

jpeg(filename="stats/Nnodes70_hist.jpeg")
hist(Nnode70_perc, main = "percentage of nodes BS>=70")
dev.off()

jpeg(filename="stats/Nnodes_vs_Nnodes70.jpeg")
plot(Nnode_perc, Nnode70_perc, main = "tree resolution", xlab = "% nodes resolved", ylab = "%nodes w/ BS>=70")
dev.off()

col = rep("black", nrow(alnstats))
col[Parsimony_informative_sites<(No_of_taxa-2)] <- "red"

pch = rep(1, nrow(alnstats))
pch[Parsimony_informative_sites<(No_of_taxa-2)] <- 16

jpeg(filename="stats/PIC_vs_Nnodes.jpeg")
plot(Parsimony_informative_sites, Nnode_perc, main ="Variation vs. resolution", xlab="No. PIC", ylab="% nodes resolved", col=col, pch=pch)
dev.off()

jpeg(filename="stats/PIC_vs_Nnodes70.jpeg")
plot(Parsimony_informative_sites, Nnode70_perc, main ="Variation vs. high resolution", xlab="No. PIC", ylab="% nodes with BS>=70", col=col, pch=pch)
dev.off()

for(name in alnstats[Parsimony_informative_sites<(No_of_taxa-2),"Alignment_name"]){
	print(name)
	tree <- read.tree(paste("raxmltrees/",name,".raxml.support", sep=""))
	jpeg(paste("stats/trees/",name,".jpg",sep=""))
	plot(tree,show.node.label=TRUE)
	dev.off()
}

write.table(alnstats, file="stats/stats.csv", sep=";", dec=",")