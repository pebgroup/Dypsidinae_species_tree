wd <- "/Users/au265104/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/gis07.st.client.au.dk – SFTP/data_vol/wolf/Dypsis"
setwd(wd)

library(ape)

# read name translation table (SECAPR No. to tip name)
rename <- read.table("rename.csv", sep=";", colClasses = "character")

# read species trees
astral_tree <- read.tree("final_tree_nofilter/astral/astral_tree.tre")
astral_tree_renamed <- read.tree("final_tree_nofilter/astral/astral_tree_renamed.tre")

# read gene trees
gtnames = list.files("final_tree_nofilter/iqtree", pattern="*.treefile", full.names=FALSE)
gts <- list()
for(i in 1:length(gtnames)){
  gts[[i]] <- read.tree(paste("final_tree_nofilter/iqtree/", gtnames[i], sep=""))
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

for(i in 1:length(gts)){
  tree = gts[[i]]
  og = c("1011","1012")
  og = og[og %in% tree$tip.label]
  if(i == 12){ og = c("1011") }
  tree = root(tree, outgroup=og)
  plot(tree)
}
  
