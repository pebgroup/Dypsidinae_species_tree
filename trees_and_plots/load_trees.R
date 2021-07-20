library(ape)

data_dir <- "/Users/au265104/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/gis07.st.client.au.dk – SFTP/data_vol/wolf/Dypsis"

# read name translation table (SECAPR No. to tip name)
rename <- read.table(paste(data_dir, "/rename.csv", sep=""), sep=";", colClasses = "character")

# read figurename translation table (SECAPR No. to figure name)
figurename <- read.table(paste(data_dir, "/figurenames2.csv", sep=""), sep=";", colClasses = "character")
figurename_idx <- figurename$V2
names(figurename_idx) <- figurename$V1

# read species trees
astral_tree <- read.tree(paste(data_dir, "/final_tree_nofilter/astral/astral_tree.tre", sep=""))
#astral_tree_renamed <- read.tree("final_tree_nofilter/astral/astral_tree_renamed.tre")
astral_tree_EN <- read.tree(paste(data_dir, "/final_tree_nofilter/astral/astral_tree_full_annot.tre", sep=""))
  
# read gene trees
gtnames = list.files(paste(data_dir, "/final_tree_nofilter/iqtree", sep=""), pattern="*.treefile", full.names=FALSE)
#gtnames = list.files("sandbox/done", pattern="*.treefile", full.names=FALSE)
gts <- list()
for(i in 1:length(gtnames)){
  gts[[i]] <- read.tree(paste(data_dir, "/final_tree_nofilter/iqtree/", gtnames[i], sep=""))
  #gts[[i]] <- read.tree(paste("sandbox/done/", gtnames[i], sep=""))
}
rm(i)