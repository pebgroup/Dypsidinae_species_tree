library(phytools)
library(geiger)
library(stringr)
library(party)

setwd("/Users/au265104/OneDrive - Aarhus Universitet/ANALYSIS/Dypsis/Dypsidinae_species_tree/trees_and_plots")

# Assuming script location within Git repo as working directory
source("functions.R")
source("load_trees.R")

traitData3 <- read.csv2("~/OneDrive - Aarhus Universitet/PROJECTS/65 Dypsis systematics paper/~manuscript/figures/trait figure/DypsidinaeTraitData+WLEv2.csv",header = TRUE)
traitData3 <- traitData3[!traitData3$SpecName=="",]
traitData3$SpecName[traitData3$SpecName == "Dypsis Leucomalla"] <- "Dypsis leucomalla"
traitData3 <- traitData3[!traitData3$SpecName=="Loxococcus rupicola",] #exclude Loxococcus

# split into species that are in the tree (this is what most of the remainder of the script works with)  and species that are NOT in the tree for adding to one of the figures
traitData3_notintree <- traitData3[traitData3$NotInTree == 1,1:45]
traitData3_complete <- traitData3 # keeping "traitData3" for the species in the tree to keep rest of the script working
traitData3 <- traitData3[traitData3$NotInTree != 1,1:45] 

# manually ladderize astral_tree_to_ladderize.tre in FigTree and save as astral_tree_lad.tre (ladderizing in R messed up the plot)
# tree <- astral_tree
# tree$edge.length[is.na(tree$edge.length)] <- 0.001
# write.tree(tree,paste(data_dir, "/final_tree_nofilter/astral/astral_tree_to_ladderize.tre", sep="")) 

tree <- read.tree(paste(data_dir, "/final_tree_nofilter/astral/astral_tree_lad_dec.tre", sep=""))
# remove outgroup and redundant samples (where more than one individual per species)
tree <- drop.tip(tree, c("0194", "0196", "0199", "0202", "0204", "1012", "1011"))

# match tree tips to trait dataset
##################################

rownames(rename) <- rename$V1

names <- rename[tree$tip.label,"V2"]
names[names == "Dypsis-heteromorpha-JD7822_S62_L001"] = "Dypsis-sp-heteromorpha-JD7822_S62_L001"
names[names == "Dypsis-mirabilis-florencei_S64_L001"] = "Dypsis-sp-mirabilis-florencei_S64_L001"
names[names == "Dypsis-onilahensis-Isalo-form_S56_L001"] = "Dypsis-sp-onilahensis-Isalo-form_S56_L001"
names[names == "Lemurophoenix-sp-Nov-SBL471"] = "Lemurophoenix-laevis-SBL471" 
names[names == "Dypsis-sp-nov-sira_S7_L001"] = "Dypsis-mijoroana_S7_L001"
names[names == "Dypsis-sp-nov-aff-lutea_S46_L001" ] = "Dypsis-aurantiaca_S46_L001" 
names <- str_replace_all(names,"-","_")
names2 <- c()
for(i in strsplit(names,"_")){
  names2 <- c(names2, paste(i[1:2], collapse="_"))
}
names2[names2 == "Dypsis_sp"] <- paste("Dypsis_sp", 1:6, sep="_")

rename2 <- figurename[1:174,]
rownames(rename2) <- rename2$V1
figurenames <- rename2[tree$tip.label,"V2"]
names(figurenames) <- names2

traitData3$figurenames <- figurenames[str_replace(traitData3$SpecName," ","_")]

# rename tree with figure names
tree$tip.label <- figurenames

# add column with new genera
head(traitData3)
novgen <- c()
for(ng in strsplit(x = traitData3$figurenames, " ")){
  novgen <- c(novgen, ng[1][1])
}
novgen[novgen == "C."] <- "Chrysalidocarpus"
novgen[novgen == "D."] <- "Dypsis"
novgen[novgen == "V."] <- "Vonitra"
novgen[novgen == "L."] <- "Lemurophoenix"
traitData3$NovGen <- novgen


# Tree with traits: Figure S1
############

#figurepath <- "/Users/au265104/OneDrive - Aarhus Universitet/PROJECTS/65 Dypsis systematics paper/~manuscript/figures/trait figure"
figurepath <- "/Users/au265104/OneDrive - Aarhus Universitet/PROJECTS/65 Dypsis systematics paper/~manuscript/figures/trait figure/new"

pdf(paste(figurepath, "/traitplot_newdata_newnames.pdf", sep=""), width=8.3, height = 11.7)

n <- length(tree$tip.label) # number of tips in the tree
plot.phylo(tree, label.offset=5.5, align.tip.label = T, cex = 0.45)

xoffset = 9.3
increment = .6
scf = 1.8

#Max stem heigth
char <- traitData3$MaxStemHeight_m
names(char) <- traitData3$figurenames
char <- char[tree$tip.label] #orders the names
points(x=rep(xoffset, n), y=1:n, pch=19, col="brown", cex=scf*char/max(char, na.rm=T))

#Max stem diam
xoffset <- xoffset + increment
char <- traitData3$MaxStemDia_cm
names(char) <- traitData3$figurenames
char <- char[tree$tip.label] #orders the names
points(x=rep(xoffset, n), y=1:n, pch=19, col="blue", cex=scf*char/max(char, na.rm=T))

#Max blade length
xoffset <- xoffset + increment
char <- traitData3$Max_Blade_Length_m
names(char) <- traitData3$figurenames
char <- char[tree$tip.label] #orders the names
points(x=rep(xoffset, n), y=1:n, pch=19, col="purple", cex=scf*char/max(char, na.rm=T))

#Fruit volume
xoffset <- xoffset + increment
char <- log(4/3 * pi * traitData3$AverageFruitLength_cm * traitData3$AverageFruitWidth_cm^2)
names(char) <- traitData3$figurenames
char <- char[tree$tip.label] #orders the names
points(x=rep(xoffset, n), y=1:n, pch=19, col="blue4", cex=scf*char/max(char, na.rm=T))

#Peduncle length
xoffset <- xoffset + increment
char <- traitData3$MaxPeduncleLength_cm
names(char) <- traitData3$figurenames
char <- char[tree$tip.label] #orders the names
points(x=rep(xoffset, n), y=1:n, pch=19, col="forestgreen", cex=scf*char/max(char, na.rm=T))

#Inflorescence length
xoffset <- xoffset + increment
char <- traitData3$MaxPeduncleLength_cm + traitData3$MaxRachisLength_cm
names(char) <- traitData3$figurenames
char <- char[tree$tip.label] #orders the names
points(x=rep(xoffset, n), y=1:n, pch=19, col="red", cex=scf*char/max(char, na.rm=T))

#Stamen number 
xoffset <- xoffset + increment
symbol <- rep(19, n)
char2 <- as.character(traitData3$StamenNumber)
names(char2) <- traitData3$figurenames
char2 <- char2[tree$tip.label]
symbol[char2=="3"] <- 1
symbol[!char2 %in% c("6", "3")] <- NA
symbol[char2=="(3) or rarely 1-2"] <- 1 
symbol[char2=="52-59"] <- 8
symbol[char2=="33-41"] <- 8
points(x=rep(xoffset, n), y=1:n, pch=symbol, cex=.6, col="orange")

#Antepetalous stamens 
xoffset <- xoffset + increment
symbol <- rep(19, n)
char2 <- as.character(traitData3$Antepetalous)
names(char2) <- traitData3$figurenames
char2 <- char2[tree$tip.label]
symbol[char2=="0"] <- 1
symbol[!char2 %in% c("0", "1")] <- NA
points(x=rep(xoffset, n), y=1:n, pch=symbol, cex=.6, col="turquoise3")

#Ruminate endosperm 
xoffset <- xoffset + increment
symbol <- rep(19, n)
char2 <- as.character(traitData3$RuminateEndosperm)
names(char2) <- traitData3$figurenames
char2 <- char2[tree$tip.label]
symbol[char2=="0"] <- 1
symbol[!char2 %in% c("0", "1")] <- NA
points(x=rep(xoffset, n), y=1:n, pch=symbol, cex=.6, col="grey31")

legend("topleft", c("Max. stem height", "Max. stem diameter", "Max. blade length", "Fruit volume", "Peduncle length", "Inflorescence length", "3 stamens", "6 stamens", ">30 stamens", "Antesepalous stamens", "Antepetalous stamens", "Smooth endosperm", "Ruminate endosperm"), 
       col = c("brown", "blue","purple","blue4", "forestgreen", "red", "orange", "orange", "orange", "turquoise3", "turquoise3", "grey31", "grey31"),
       text.col = "black", pch = c(19,19,19,19,19,19,1,19,8,1,19,1,19), bg = "white")

dev.off()


# Classification
################

classdata <- traitData3[traitData3$NovGen %in% c("Dypsis", "Chrysalidocarpus"),]
classdata$NovGen <- as.factor(classdata$NovGen)

classdata <- classdata[,c("NovGen", "MaxStemHeight_m", "MaxStemDia_cm", "MinStemDia_cm", "MaxLeafNumber", "Max_Blade_Length_m", "Max_Petiole_length_m", "AverageFruitLength_cm", "AverageFruitWidth_cm"
                        , "MaxPeduncleLength_cm", "MaxRachisLength_cm")]

#combine rachis and lamina length into leaf length
classdata$MaxLeafLength <- classdata$Max_Blade_Length_m + classdata$Max_Petiole_length_m
classdata <- classdata[,!colnames(classdata) %in% c("Max_Blade_Length_m", "Max_Petiole_length_m")]

# Classification (Figure S2)
d_vs_c <- ctree(NovGen ~ MaxStemHeight_m + MaxStemDia_cm + MaxLeafNumber + MaxLeafLength + AverageFruitLength_cm + AverageFruitWidth_cm + MaxPeduncleLength_cm + MaxRachisLength_cm,  data = classdata)
plot(d_vs_c)

# Minimum stem height
d_vs_c <- ctree(NovGen ~ MaxStemHeight_m + MinStemDia_cm + MaxLeafNumber + MaxLeafLength + AverageFruitLength_cm + AverageFruitWidth_cm + MaxPeduncleLength_cm + MaxRachisLength_cm,  data = classdata)
plot(d_vs_c)


# Figure 6
##########

figurepath <- "/Users/au265104/OneDrive - Aarhus Universitet/PROJECTS/65 Dypsis systematics paper/~manuscript/figures/classification figs"

plotdata <- traitData3#[traitData3$NovGen %in% c("Dypsis", "Chrysalidocarpus"),]
plotdata$NovGen <- as.factor(plotdata$NovGen)
#plotdata[plotdata$figurenames == "D. makirae", "Max_Blade_Length_m"] <- 0.15 + 0.25 # Petiole + Rachis, Rakotoarinivo et al. 2009 PALMS

plotdata <- plotdata[,c("NovGen", "MaxStemHeight_m", "MaxStemDia_cm", "MaxLeafNumber", "Max_Blade_Length_m", "Max_Petiole_length_m", "AverageFruitLength_cm", "AverageFruitWidth_cm"
                        , "MaxPeduncleLength_cm", "MaxRachisLength_cm")]

#combine rachis and lamina length into leaf length
plotdata$MaxLeafLength <- plotdata$Max_Blade_Length_m + plotdata$Max_Petiole_length_m
plotdata <- plotdata[,!colnames(plotdata) %in% c("Max_Blade_Length_m", "Max_Petiole_length_m")]

# exclude NA species
plotdata <- plotdata[!(is.na(plotdata$MaxLeafLength) | is.na(plotdata$MaxStemDia_cm)),]

colorpal <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
pchpal <- c(16,2,3,4,5,6)

pdf(paste(figurepath, "/biplot_leaflength_all_new.pdf", sep=""), width=1.5*6.93, height = 1.5*6.93/2)
#pdf(paste(figurepath, "/biplot_leaflength_all_txt.pdf", sep=""), width=2*1.5*6.93, height = 2*1.5*6.93/2)
par(mfrow=c(1,2))

plotdata_zoom <- plotdata[plotdata$MaxStemDia_cm<10 & plotdata$MaxLeafLength<2.5,]

col = rep(colorpal[1], nrow(plotdata))
col[plotdata$NovGen == "Dypsis"] <- colorpal[2]
col[plotdata$NovGen == "Vonitra"] <- colorpal[3]
col[plotdata$NovGen == "Marojejya"] <- colorpal[4]
col[plotdata$NovGen == "Masoala"] <- colorpal[5]
col[plotdata$NovGen == "Lemurophoenix"] <- colorpal[6]

pch = rep(pchpal[1], nrow(plotdata))
pch[plotdata$NovGen == "Dypsis"] <- pchpal[2]
pch[plotdata$NovGen == "Vonitra"] <- pchpal[3]
pch[plotdata$NovGen == "Marojejya"] <- pchpal[4]
pch[plotdata$NovGen == "Masoala"] <- pchpal[5]
pch[plotdata$NovGen == "Lemurophoenix"] <- pchpal[6]
pch[plotdata$MaxStemDia_cm == 2 & plotdata$MaxLeafLength == 1.2 + 0.35] <- 1

plot(plotdata$MaxLeafLength, plotdata$MaxStemDia_cm, col = col, cex=1, xlab="Maximum leaf length [m]", ylab="Maximum stem diameter [cm]", pch=pch)
abline(v=1.65)
abline(h=6)
#text(plotdata$MaxLeafLength, plotdata$MaxStemDia_cm, str_replace(plotdata$SpecName, "Dypsis ", ""), cex=0.4)
#text(traitData3_notintree$Max_Blade_Length_m + traitData3_notintree$Max_Petiole_length_m, traitData3_notintree$MaxStemDia_cm, str_replace(traitData3_notintree$SpecName, "Dypsis ", ""), cex = 0.4, col = "red")

lines(c(-0.1,max(plotdata_zoom$MaxLeafLength, na.rm=T)), c(max(plotdata_zoom$MaxStemDia_cm, na.rm=T),max(plotdata_zoom$MaxStemDia_cm, na.rm=T)), lty=3)
lines(c(max(plotdata_zoom$MaxLeafLength, na.rm=T),max(plotdata_zoom$MaxLeafLength, na.rm=T)), c(max(plotdata_zoom$MaxStemDia_cm, na.rm=T),-2), lty=3)
legend("topright", legend=c("Chrysalidocarpus", "Dypsis", "Vonitra", "Marojejya", "Masoala", "Lemurophoenix"), pch=pchpal, col=colorpal[1:6], cex=.75, text.font = 3) #, bty="n"

col = rep(colorpal[1], nrow(plotdata_zoom))
col[plotdata_zoom$NovGen == "Dypsis"] <- colorpal[2]
col[plotdata_zoom$NovGen == "Vonitra"] <- colorpal[3]

pch = rep(pchpal[1], nrow(plotdata_zoom))
pch[plotdata_zoom$NovGen == "Dypsis"] <- pchpal[2]
pch[plotdata_zoom$NovGen == "Vonitra"] <- pchpal[3]
pch[plotdata_zoom$MaxStemDia_cm == 2 & plotdata_zoom$MaxLeafLength == 1.2 + 0.35] <- 1


plot(plotdata_zoom$MaxLeafLength, plotdata_zoom$MaxStemDia_cm, col = col, cex=1, xlab="Maximum leaf length [m]", ylab="Maximum stem diameter [cm]", pch=pch)
abline(v=1.65)
abline(h=6)
#text(plotdata_zoom$MaxLeafLength, plotdata_zoom$MaxStemDia_cm, str_replace(plotdata_zoom$SpecName, "Dypsis ", ""), cex=0.4)
#text(traitData3_notintree$Max_Blade_Length_m + traitData3_notintree$Max_Petiole_length_m, traitData3_notintree$MaxStemDia_cm, str_replace(traitData3_notintree$SpecName, "Dypsis ", ""), cex = 0.4, col = "red")
dev.off()


# Leaf length and MINIMUM stem diameter with all genera
##########################

plotdata <- traitData3#[traitData3$NovGen %in% c("Dypsis", "Chrysalidocarpus"),]
plotdata$NovGen <- as.factor(plotdata$NovGen)
#plotdata[plotdata$figurenames == "D. makirae", "Max_Blade_Length_m"] <- 0.15 + 0.25 # Petiole + Rachis, Rakotoarinivo et al. 2009 PALMS

plotdata <- plotdata[,c("SpecName", "NovGen", "MaxStemHeight_m", "MinStemDia_cm", "MaxLeafNumber", "Max_Blade_Length_m", "Max_Petiole_length_m",  "AverageFruitLength_cm", "AverageFruitWidth_cm"
                        , "MaxPeduncleLength_cm", "MaxRachisLength_cm", "StamenNumber", "InflorescenceBranchingOrder_norm")]

#combine rachis and lamina length into leaf length
plotdata$MaxLeafLength <- plotdata$Max_Blade_Length_m + plotdata$Max_Petiole_length_m
plotdata <- plotdata[,!colnames(plotdata) %in% c("Max_Blade_Length_m", "Max_Petiole_length_m")]

plotdata[is.na(plotdata[,"InflorescenceBranchingOrder_norm"]),"InflorescenceBranchingOrder_norm"] <- -99

# exclude NA species
plotdata <- plotdata[!(is.na(plotdata$MaxLeafLength) | is.na(plotdata$MinStemDia_cm)),]

#classdata <- plotdata[plotdata$NovGen %in% c("Dypsis", "Chrysalidocarpus"),2:11]
#classdata$NovGen <- droplevels(classdata$NovGen)
#d_vs_c <- ctree(NovGen ~ .,  data = classdata)
#plot(d_vs_c)

colorpal <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
pchpal <- c(16,2,3,4,5,6)

#pdf(paste(figurepath, "/biplot_leaflength_all_mindiam.pdf", sep=""), width=1.5*6.93, height = 1.5*6.93/2)
pdf(paste(figurepath, "/biplot_leaflength_all_mindiam_txt.pdf", sep=""), width=2*1.5*6.93, height = 2*1.5*6.93/2)
par(mfrow=c(1,2))

plotdata_zoom <- plotdata[plotdata$MinStemDia_cm<10 & plotdata$MaxLeafLength<2.5,]

col = rep(colorpal[1], nrow(plotdata))
col[plotdata$NovGen == "Dypsis"] <- colorpal[2]
col[plotdata$NovGen == "Vonitra"] <- colorpal[3]
col[plotdata$NovGen == "Marojejya"] <- colorpal[4]
col[plotdata$NovGen == "Masoala"] <- colorpal[5]
col[plotdata$NovGen == "Lemurophoenix"] <- colorpal[6]

pch = rep(pchpal[1], nrow(plotdata))
pch[plotdata$NovGen == "Dypsis"] <- pchpal[2]
pch[plotdata$NovGen == "Vonitra"] <- pchpal[3]
pch[plotdata$NovGen == "Marojejya"] <- pchpal[4]
pch[plotdata$NovGen == "Masoala"] <- pchpal[5]
pch[plotdata$NovGen == "Lemurophoenix"] <- pchpal[6]
pch[plotdata$MinStemDia_cm == 2 & plotdata$MaxLeafLength == 1.2 + 0.35] <- 1

plot(plotdata$MaxLeafLength, plotdata$MinStemDia_cm, col = col, cex=1, xlab="Maximum leaf length [m]", ylab="Minimum stem diameter [cm]", pch=pch)
abline(v=1.75)
abline(h=2)
lines(c(-0.1,max(plotdata_zoom$MaxLeafLength, na.rm=T)), c(max(plotdata_zoom$MinStemDia_cm, na.rm=T),max(plotdata_zoom$MinStemDia_cm, na.rm=T)), lty=3)
lines(c(max(plotdata_zoom$MaxLeafLength, na.rm=T),max(plotdata_zoom$MaxLeafLength, na.rm=T)), c(max(plotdata_zoom$MinStemDia_cm, na.rm=T),-2), lty=3)
txtcol <- rep("black", nrow(plotdata))
txtcol[plotdata$StamenNumber == 3 | plotdata$InflorescenceBranchingOrder_norm == 0] <- "blue"
#text(plotdata$MaxLeafLength, plotdata$MinStemDia_cm, str_replace(plotdata$SpecName, "Dypsis ", ""), cex=0.4, col= txtcol)
#text(traitData3_notintree$Max_Blade_Length_m + traitData3_notintree$Max_Petiole_length_m, traitData3_notintree$MinStemDia_cm, str_replace(traitData3_notintree$SpecName, "Dypsis ", ""), cex = 0.4, col = "red")
legend("topright", legend=c("Chrysalidocarpus", "Dypsis", "Vonitra", "Marojejya", "Masoala", "Lemurophoenix"), pch=pchpal, col=colorpal[1:6], cex=.75, text.font = 3) #, bty="n"

col = rep(colorpal[1], nrow(plotdata_zoom))
col[plotdata_zoom$NovGen == "Dypsis"] <- colorpal[2]
col[plotdata_zoom$NovGen == "Vonitra"] <- colorpal[3]

pch = rep(pchpal[1], nrow(plotdata_zoom))
pch[plotdata_zoom$NovGen == "Dypsis"] <- pchpal[2]
pch[plotdata_zoom$NovGen == "Vonitra"] <- pchpal[3]
pch[plotdata_zoom$MinStemDia_cm == 2 & plotdata_zoom$MaxLeafLength == 1.2 + 0.35] <- 1


plot(plotdata_zoom$MaxLeafLength, plotdata_zoom$MinStemDia_cm, col = col, cex=1, xlab="Maximum leaf length [m]", ylab="Minimum stem diameter [cm]", pch=pch)
abline(v=1.75)
abline(h=2)
txtcol <- rep("black", nrow(plotdata_zoom))
txtcol[plotdata_zoom$StamenNumber == 3 | plotdata_zoom$InflorescenceBranchingOrder_norm == 0] <- "blue"
#text(plotdata_zoom$MaxLeafLength, plotdata_zoom$MinStemDia_cm, str_replace(plotdata_zoom$SpecName, "Dypsis ", ""), cex=0.4, col=txtcol)
#text(traitData3_notintree$Max_Blade_Length_m + traitData3_notintree$Max_Petiole_length_m, traitData3_notintree$MinStemDia_cm, str_replace(traitData3_notintree$SpecName, "Dypsis ", ""), cex = 0.4, col = "red")

dev.off()

traitData3[traitData3$NovGen=="Dypsis",] -> dypsis
traitData3[traitData3$NovGen=="Chrysalidocarpus",] -> chrysa
hist(dypsis$MaxStemDia_cm, breaks=0:31, xlab="Maximum stem diameter", main="Dypsis s.str.")
sum(dypsis$MaxStemDia_cm <= 2, na.rm=T)
nrow(dypsis)-sum(is.na(dypsis$MaxStemDia_cm))

hist(chrysa$MinStemDia_cm, breaks=0:45)

