#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(ape)
library(phytools)
library(phangorn)
library(RColorBrewer)

spptree = read.tree("astral_tree.tre")

gtnames = list.files("../length_filter/iqtree", pattern="*.treefile", full.names=FALSE)

# Define groups for which to check monophyly

dypsis = spptree$tip.label[Descendants(spptree, findMRCA(spptree, c("0060", "0035")), "tips")[[1]]]

chrysalido = spptree$tip.label[Descendants(spptree, findMRCA(spptree, c("0078", "0174")), "tips")[[1]]]

maro = spptree$tip.label[Descendants(spptree, findMRCA(spptree, c("0198", "0205")), "tips")[[1]]]

vonitra = spptree$tip.label[Descendants(spptree, findMRCA(spptree, c("0055", "0116")), "tips")[[1]]] #0031
vonitra_large = spptree$tip.label[Descendants(spptree, findMRCA(spptree, c("0031", "0116")), "tips")[[1]]] #0031

lemuro = c("0191", "0185")

maso = c("0077", "0076", "0187")

taxa = c("maso", "lemuro", "maro", "dypsis", "chrysalido", "vonitra")

labels = c("Masoala", "Lemurophoenix", "Marojejya", "Dypsis s.str.", "Chrysalidocarpus", "Vonitra")
names(labels) = taxa

problems = rep(0,length(gtnames)*7)
dim(problems) = c(length(gtnames),7)
problems = as.data.frame(problems)
colnames(problems) <- c(taxa, "vonitra_large")
rownames(problems) <- gtnames

# plotting the distribution of support

jpeg(filename="Monophyly_support.jpg", width = 800, height = 1200, quality=100, res=200)
par(mfrow=c(3,2))

for(t in taxa){

	sup = c() # vector of support values in trees where the taxon is monophyletic
	nomo = 0 # number of trees in which testing is not possible due to taxon sampling
	
	taxon = get(t)

	for(i in 1:length(gtnames)){
		tree = read.tree(paste("../length_filter/iqtree/", gtnames[i], sep=""))
		og = c("1011","1012")
		og = og[og %in% tree$tip.label]
		tree = root(tree, outgroup=og)

		clade = taxon

		clade = clade[clade %in% tree$tip.label]

		if(length(clade) > 1){ #if monophyly-check even possible
		  if(length(Descendants(tree, findMRCA(tree, clade), "tips")[[1]]) == length(clade)){ # is monophyletic
			sup = c(sup, tree$node.label[findMRCA(tree, clade)-length(tree$tip.label)])
		  } else {
		  	  problems[gtnames[i],t] <- 1
			  #colour = rep("black", length(tree$tip.label))
			  #colour[tree$tip.label %in% clade] = "red"
			  #jpeg(filename=paste("dypsis/",gtnames[i],".jpeg",sep=""), width=800, height=1400)
			  #plot(tree, tip.color = colour, main=gtnames[i])
			  #dev.off()
		  }
		} else {nomo = nomo+1}
	}

	sup = as.numeric(sup)

	# calculate for each bootstrap threshold the number of genetrees in which the clade is supported above that threshold
	sup_freq = c()

	for(j in 30:100){
		sup_freq <- c(sup_freq, sum(sup >= j))
	}


	plot(1, type="n", xlim=c(30,100), ylim=c(1,length(gtnames)+10), ylab="number of gene trees", xlab="BS", main=labels[t])
	lines(30:100, sup_freq, lwd=2)
	abline(h=length(gtnames)-nomo, lty=3)
	
	if(t == "vonitra"){
	
		sup = c() 
		nomo = 0 
	
		taxon = get("vonitra_large")

		for(i in 1:length(gtnames)){
			tree = read.tree(paste("../length_filter/iqtree/", gtnames[i], sep=""))
			og = c("1011","1012")
			og = og[og %in% tree$tip.label]
			tree = root(tree, outgroup=og)

			clade = taxon

			clade = clade[clade %in% tree$tip.label]

			if(length(clade) > 1){ #if monophyly-check even possible
			  if(length(Descendants(tree, findMRCA(tree, clade), "tips")[[1]]) == length(clade)){ # is monophyletic
				sup = c(sup, tree$node.label[findMRCA(tree, clade)-length(tree$tip.label)])
			  } else {
			      problems[gtnames[i],"vonitra_large"] <- 1
				  #colour = rep("black", length(tree$tip.label))
				  #colour[tree$tip.label %in% clade] = "red"
				  #jpeg(filename=paste("dypsis/",gtnames[i],".jpeg",sep=""), width=800, height=1400)
				  #plot(tree, tip.color = colour, main=gtnames[i])
				  #dev.off()
			  }
			} else {nomo = nomo+1}
		}

		sup = as.numeric(sup)

		# calculate for each bootstrap threshold the number of genetrees in which the clade is supported above that threshold
		sup_freq = c()

		for(j in 30:100){
			sup_freq <- c(sup_freq, sum(sup >= j))
		}

		lines(30:100, sup_freq, lwd=2, col="grey")
		abline(h=length(gtnames)-nomo, lty=3, col="grey")
	
	}
	

}
dev.off()

write.csv(problems,"problematic_genetrees.csv")

# for the three possible resolutions of Dypsis-Marojejya

res1 = c(dypsis, chrysalido)
res2 = c(dypsis, maro)
res3 = c(maro, chrysalido)

taxa = c("res1", "res2", "res3")

labels = c("Dypsis s.str. + Chrysalidocarpus", "Dypsis s.str. + Marojejya", "Chrysalidocarpus + Marojejya")
names(labels) = taxa

cols = brewer.pal(3, "Set3")
names(cols) = taxa

jpeg(filename="Alternative_resolutions.jpg", width = 800, height = 800, quality=100, res=200)
plot(1, type="n", xlim=c(30,100), ylim=c(1,length(gtnames)+10), ylab="number of gene trees", xlab="BS", main="Alternative resolutions")

sup_old = rep(0, 71)
for(t in taxa){

	sup = c() # vector of support values in trees where the taxon is monophyletic
	nomo = 0 # number of trees in which testing is not possible due to taxon sampling
	
	taxon = get(t)

	for(i in 1:length(gtnames)){
		tree = read.tree(paste("../length_filter/iqtree/", gtnames[i], sep=""))
		og = c("1011","1012")
		og = og[og %in% tree$tip.label]
		tree = root(tree, outgroup=og)

		clade = taxon

		clade = clade[clade %in% tree$tip.label]

		if(length(clade) > 1){ #if monophyly-check even possible
		  if(length(Descendants(tree, findMRCA(tree, clade), "tips")[[1]]) == length(clade)){ # is monophyletic
			sup = c(sup, tree$node.label[findMRCA(tree, clade)-length(tree$tip.label)])
		  } else {
			  #colour = rep("black", length(tree$tip.label))
			  #colour[tree$tip.label %in% clade] = "red"
			  #jpeg(filename=paste("dypsis/",gtnames[i],".jpeg",sep=""), width=800, height=1400)
			  #plot(tree, tip.color = colour, main=gtnames[i])
			  #dev.off()
		  }
		} else {nomo = nomo+1}
	}

	sup = as.numeric(sup)

	# calculate for each bootstrap threshold the number of genetrees in which the clade is supported above that threshold
	sup_freq = c()

	for(j in 30:100){
		sup_freq <- c(sup_freq, sum(sup >= j))
	}
	

	x = c(30:100, 100:30, 30)
	y = c(sup_old, rev(sup_freq+sup_old), sup_freq[1]+sup_old[1])

	#print(x)
	#print(y)

	#lines(30:100, sup_freq, lwd=2)
	polygon(x,y,col=cols[t], border=NA)
	abline(h=length(gtnames)-nomo, lty=3)
	sup_old = sup_old + sup_freq
	
}

text(32,22,labels[1], adj=0, cex=.8)
text(32,57,labels[2], adj=0, cex=.8)
text(32,95,labels[3], adj=0, cex=.8)
dev.off()