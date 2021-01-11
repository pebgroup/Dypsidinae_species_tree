library(ape)

# read name translation table (SECAPR No. to tip name)
rename <- read.table("rename.csv", sep=";", colClasses = "character")

# read figurename translation table (SECAPR No. to figure name)
figurename <- read.table("figurenames.csv", sep=";", colClasses = "character")
figurename_idx <- figurename$V2
names(figurename_idx) <- figurename$V1

# read species trees
astral_tree <- read.tree("final_tree_nofilter/astral/astral_tree.tre")
astral_tree_renamed <- read.tree("final_tree_nofilter/astral/astral_tree_renamed.tre")

# read gene trees
gtnames = list.files("final_tree_nofilter/iqtree", pattern="*.treefile", full.names=FALSE)
#gtnames = list.files("sandbox/done", pattern="*.treefile", full.names=FALSE)
gts <- list()
for(i in 1:length(gtnames)){
  gts[[i]] <- read.tree(paste("final_tree_nofilter/iqtree/", gtnames[i], sep=""))
  #gts[[i]] <- read.tree(paste("sandbox/done/", gtnames[i], sep=""))
}
rm(i)