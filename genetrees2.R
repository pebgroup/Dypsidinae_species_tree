wd <- "/Users/au265104/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/gis07.st.client.au.dk – SFTP/data_vol/wolf/Dypsis"
setwd(wd)

library(ape)
library(phytools)
library(phangorn)
library(RColorBrewer)

source("/Users/au265104/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/gis07.st.client.au.dk – SFTP/home/au265104/scripts/dypsidinae/load_trees.R")

#load_trees.R
#
# # read name translation table (SECAPR No. to tip name)
# rename <- read.table("rename.csv", sep=";", colClasses = "character")
# 
# # read figurename translation table (SECAPR No. to figure name)
# figurename <- read.table("figurenames.csv", sep=";", colClasses = "character")
# figurename_idx <- figurename$V2
# names(figurename_idx) <- figurename$V1
# 
# # read species trees
# astral_tree <- read.tree("final_tree_nofilter/astral/astral_tree.tre")
# astral_tree_renamed <- read.tree("final_tree_nofilter/astral/astral_tree_renamed.tre")
# 
# # read gene trees
# gtnames = list.files("final_tree_nofilter/iqtree", pattern="*.treefile", full.names=FALSE)
# #gtnames = list.files("sandbox/done", pattern="*.treefile", full.names=FALSE)
# gts <- list()
# for(i in 1:length(gtnames)){
#   gts[[i]] <- read.tree(paste("final_tree_nofilter/iqtree/", gtnames[i], sep=""))
#   #gts[[i]] <- read.tree(paste("sandbox/done/", gtnames[i], sep=""))
# }
# rm(i)

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

###################
# Check monophyly #
###################

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

problems = rep(0,length(gtnames)*7)
dim(problems) = c(length(gtnames),7)
problems = as.data.frame(problems)
colnames(problems) <- c(taxa, "vonitra_large")
rownames(problems) <- gtnames

# plotting the distribution of support

jpeg(filename="/Users/au265104/Desktop/Monophyly_support.jpg", width = 800, height = 1200, quality=100, res=200)
par(mfrow=c(3,2))

for(t in taxa){
  
  sup = c() # vector of support values in trees where the taxon is monophyletic
  nomo = 0 # number of trees in which testing is not possible due to taxon sampling
  
  taxon = get(t)
  
  for(i in 1:length(gts)){
    tree = gts[[i]]
    og = c("1011","1012")
    og = og[og %in% tree$tip.label]
    if(i == 12){ og = c("1011") }
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
    
    for(i in 1:length(gts)){
      tree = gts[[i]]
      og = c("1011","1012")
      og = og[og %in% tree$tip.label]
      if(i == 12){ og = c("1011") }
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

#write.csv(problems,"problematic_genetrees.csv")

##################
# Plot genetrees #
##################

flagged <- row.names(problems)[rowSums(problems[,1:6])>0]

colpal = c(brewer.pal(5, "Set1"), "maroon1")

blacklist = c("0075", "0076", "0157", "0197", "0159", "0164", "2013", "2016", "0119")

titlecol = "black"
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
  #if("2051" %in% tree$tip.label) titlecol = "red" 
  if(gtnames[i] %in% flagged) titlecol = "red" 
  tree$tip.label = figurename_idx[tree$tip.label]
  plot(tree, cex=.5, tip.color = tipcol, show.node.label = TRUE)
  title(main= strsplit(gtnames[i],"_")[[1]][2], col.main = titlecol)
  titlecol = "black"
  dev.off()
}
