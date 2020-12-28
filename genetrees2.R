wd <- "/Users/au265104/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/gis07.st.client.au.dk – SFTP/data_vol/wolf/Dypsis"
setwd(wd)

library(ape)
library(phytools)
library(phangorn)
library(RColorBrewer)

# read name translation table (SECAPR No. to tip name)
rename <- read.table("rename.csv", sep=";", colClasses = "character")

# read species trees
astral_tree <- read.tree("final_tree_nofilter_v1/astral/astral_tree.tre")
astral_tree_renamed <- read.tree("final_tree_nofilter_v1/astral/astral_tree_renamed.tre")

# read gene trees
gtnames = list.files("final_tree_nofilter_v1/iqtree", pattern="*.treefile", full.names=FALSE)
gts <- list()
for(i in 1:length(gtnames)){
  gts[[i]] <- read.tree(paste("final_tree_nofilter_v1/iqtree/", gtnames[i], sep=""))
}
rm(i)

# record representation of samples across gene trees
gt_incidence <- rep(0, nrow(rename))
names(gt_incidence) <- rename$V1
for(i in 1:length(gts)){
  for(n in names(gt_incidence)){
    if(n %in% gts[[i]]$tip.label){
      gt_incidence[n] = gt_incidence[n] + 1
    }
  }
}

write.table(gt_incidence, "/Users/au265104/Desktop/gt_incidence.csv",sep="|")


# Plot genetrees

# Define groups for which to check monophyly

dypsis = astral_tree$tip.label[Descendants(astral_tree, findMRCA(astral_tree, c("0060", "0035")), "tips")[[1]]]

chrysalido = astral_tree$tip.label[Descendants(astral_tree, findMRCA(astral_tree, c("0078", "0174")), "tips")[[1]]]

maro = astral_tree$tip.label[Descendants(astral_tree, findMRCA(astral_tree, c("0198", "0205")), "tips")[[1]]]

vonitra = astral_tree$tip.label[Descendants(astral_tree, findMRCA(astral_tree, c("0055", "0116")), "tips")[[1]]] #0031

vonitra_large = astral_tree$tip.label[Descendants(astral_tree, findMRCA(astral_tree, c("0031", "0116")), "tips")[[1]]] #0031

lemuro = c("0191", "0185")

maso = c("0077", "0187")

taxa = c("maso", "lemuro", "maro", "dypsis", "chrysalido", "vonitra")

labels = c("Masoala", "Lemurophoenix", "Marojejya", "Dypsis s.str.", "Chrysalidocarpus", "Vonitra")
names(labels) = taxa

colpal = c(brewer.pal(5, "Set1"), "maroon1")

blacklist = c("0075", "0076", "0157", "0197", "0159", "0164", "2013", "2016", "0119")

for(i in 1:length(gts)){
  pdf(paste("/Users/au265104/Desktop/genetrees/",gtnames[i],".pdf",sep=""), height=11.75, width=8.25)
  tree = gts[[i]]
  og = c("1011","1012")
  og = og[og %in% tree$tip.label]
  if(i == 12){ og = c("1011") }
  tree = root(tree, outgroup=og)
  # assign colours to tip labels
  tipcol = rep("black", length(tree$tip.label))
  for(j in 1:length(taxa)){
    tax = get(taxa[j])
    tipcol[tree$tip.label %in% tax] <- colpal[j]
    tipcol[tree$tip.label %in% blacklist] <- "lightgrey"
  }
  plot(tree, cex=.5, tip.color = tipcol)
  dev.off()
}

# put the gaff back in 
