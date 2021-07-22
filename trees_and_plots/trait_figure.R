library(phytools)
library(geiger)
library(stringr)

# Assuming script location within Git repo as working directory
source("functions.R")
source("load_trees.R")

traitData3 <- read.csv2("~/OneDrive - Aarhus Universitet/PROJECTS/65 Dypsis systematics paper/~manuscript/figures/trait figure/DypsidinaeTraitData+WLE.csv",header = TRUE)
traitData3 <- traitData3[!traitData3$SpecName=="",]
traitData3$SpecName[traitData3$SpecName == "Dypsis Leucomalla"] <- "Dypsis leucomalla"
traitData3 <- traitData3[!traitData3$SpecName=="Loxococcus rupicola",] #exclude Loxococcus

# manually ladderize astral_tree_to_ladderize.tre in FigTree and save as astral_tree_lad.tre (ladderizing in R messed up the plot)
# tree <- astral_tree
# tree$edge.length[is.na(tree$edge.length)] <- 0.001
# write.tree(tree,paste(data_dir, "/final_tree_nofilter/astral/astral_tree_to_ladderize.tre", sep="")) 

tree <- read.tree(paste(data_dir, "/final_tree_nofilter/astral/astral_tree_lad.tre", sep=""))
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
###############################
tree$tip.label <- figurenames

# Build plot
############

figurepath <- "/Users/au265104/OneDrive - Aarhus Universitet/PROJECTS/65 Dypsis systematics paper/~manuscript/figures/trait figure"

pdf(paste(figurepath, "/traitplot_newdata.pdf", sep=""), width=8.3, height = 11.7)

n <- length(tree$tip.label) # number of tips in the tree
plot.phylo(tree, label.offset=5.5, align.tip.label = T, cex = 0.45)

xoffset = 9
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

# Phylogenetic PCA
##################

# Prepare genera for plotting in different colours
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

groups <- as.data.frame(traitData3[,"NovGen"])
rownames(groups) <- traitData3$figurenames

# prepare subset of data for PCA
sizes <- traitData3[,c("MaxStemHeight_m", "MaxStemDia_cm", "Max_Blade_Length_m")]
rownames(sizes) <- traitData3$figurenames
sizes <- na.omit(sizes)
sizes[sizes$MaxStemHeight_m == 0,"MaxStemHeight_m"] <- 0.001

#The data should be log transformed as the data in sizes have different orders of magnitude:
sizes$MaxStemHeight_m <- log10(sizes$MaxStemHeight_m)
sizes$MaxStemDia_cm <- log10(sizes$MaxStemDia_cm)
sizes$Max_Blade_Length_m <- log10(sizes$Max_Blade_Length_m)

pca_tree <- drop.tip(tree,name.check(tree,sizes)$tree_not_data)

#Run phylo PCA 
dypsis.pca <- phyl.pca(pca_tree,sizes)

#Make PhyloPCA plot:
pdf(paste(figurepath, "/pcaplot_newdata.pdf", sep=""), width=5.1, height = 5.1)
phylomorphospace(pca_tree,scores(dypsis.pca)[,1:2],node.size=c(0,0.3),bty="n",ftype="off")
eco<-setNames(groups[,1],rownames(groups))
eco<-as.factor(eco)
ECO<-to.matrix(eco,levels(eco))
tiplabels(pie=ECO[pca_tree$tip.label,],cex=0.4,piecol=c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")[2:7])
legend("topleft",legend=levels(eco),cex=0.8,pch=21,
       pt.bg=c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")[2:7],pt.cex=1.5,bty="n")
dev.off()

#plot phyloPCA with species names (for transferring selected names to final plot):
phylomorphospace(pca_tree,scores(dypsis.pca)[,1:2])


#Trait ranges
#############

novgen <- c()
for(ng in strsplit(x = traitData3$figurenames, " ")){
  novgen <- c(novgen, ng[1][1])
}
novgen[novgen == "C."] <- "Chrysalidocarpus"
novgen[novgen == "D."] <- "Dypsis"
novgen[novgen == "V."] <- "Vonitra"
novgen[novgen == "L."] <- "Lemurophoenix"
traitData3$NovGen <- novgen

traitData3$frtvlm <- round(4/3 * pi * traitData3$AverageFruitLength_cm * traitData3$AverageFruitWidth_cm^2, 1)
traitData3$infllen <- traitData3$MaxPeduncleLength_cm + traitData3$MaxRachisLength_cm

traitData3$MaxPeduncleLength_cm

traitData3[traitData3$NovGen == "Masoala","infllen"]

traitData3[traitData3$NovGen == "Lemurophoenix","infllen"]

plot(traitData3[traitData3$NovGen == "Vonitra","infllen"])
traitData3[traitData3$NovGen == "Vonitra" & !is.na(traitData3$infllen),c("SpecName","infllen")][sort(traitData3[traitData3$NovGen == "Vonitra" & !is.na(traitData3$infllen),"infllen"], index.return=T)$ix ,]

traitData3[traitData3$NovGen == "Marojejya","infllen"]

plot(traitData3[traitData3$NovGen == "Chrysalidocarpus" & !is.na(traitData3$infllen),"infllen"])
traitData3[traitData3$NovGen == "Chrysalidocarpus" & !is.na(traitData3$infllen),c("SpecName","infllen")][sort(traitData3[traitData3$NovGen == "Chrysalidocarpus" & !is.na(traitData3$infllen),"infllen"], index.return=T)$ix ,]

quantile(traitData3[traitData3$NovGen == "Dypsis" & !is.na(traitData3$infllen),"infllen"],0.95)
plot(traitData3[traitData3$NovGen == "Dypsis" & !is.na(traitData3$infllen),"infllen"])
traitData3[traitData3$NovGen == "Dypsis" & !is.na(traitData3$infllen),c("SpecName","infllen")][sort(traitData3[traitData3$NovGen == "Dypsis" & !is.na(traitData3$infllen),"infllen"], index.return=T)$ix ,]

# Classification
################

#install.packages("party")
library(party)

plotdata <- traitData3[traitData3$NovGen %in% c("Dypsis", "Chrysalidocarpus"),]
plotdata$NovGen <- as.factor(plotdata$NovGen)
#plotdata[plotdata$figurenames == "D. makirae", "Max_Blade_Length_m"] <- 0.15 + 0.25 # Petiole + Rachis, Rakotoarinivo et al. 2009 PALMS

plotdata[is.na(plotdata$Max_Blade_Length_m),]
plotdata <- plotdata[,c("NovGen", "MaxStemHeight_m", "MaxStemDia_cm", "MaxLeafNumber", "Max_Blade_Length_m", "Max_Petiole_length_m", "AverageFruitLength_cm", "AverageFruitWidth_cm"
 , "MaxPeduncleLength_cm", "MaxRachisLength_cm")]

figurepath <- "/Users/au265104/OneDrive - Aarhus Universitet/PROJECTS/65 Dypsis systematics paper/~manuscript/figures/classification figs"

plot(plotdata[,2:ncol(plotdata)])
cor(plotdata[,2:ncol(plotdata)], use="co")

d_vs_c <- ctree(NovGen ~ .,  data = plotdata)
plot(d_vs_c)

pdf(paste(figurepath, "/biplot_newdata.pdf", sep=""), width=1.5*6.93, height = 1.5*6.93/2)
par(mfrow=c(1,2))

plotdata_zoom <- plotdata[plotdata$MaxStemDia_cm<10 & plotdata$Max_Blade_Length_m<2,]

col = rep("black", nrow(plotdata))
col[plotdata$NovGen == "Dypsis"] <- "red"
plot(plotdata$Max_Blade_Length_m, plotdata$MaxStemDia_cm, col = col, cex=1, xlab="Maximum lamina length [m]", ylab="Maximum stem diameter [cm]")
abline(v=1.30)
abline(h=6)
lines(c(-0.1,max(plotdata_zoom$Max_Blade_Length_m, na.rm=T)), c(max(plotdata_zoom$MaxStemDia_cm, na.rm=T),max(plotdata_zoom$MaxStemDia_cm, na.rm=T)), lty=3)
lines(c(max(plotdata_zoom$Max_Blade_Length_m, na.rm=T),max(plotdata_zoom$Max_Blade_Length_m, na.rm=T)), c(max(plotdata_zoom$MaxStemDia_cm, na.rm=T),-2), lty=3)

col = rep("black", nrow(plotdata_zoom))
col[plotdata_zoom$NovGen == "Dypsis"] <- "red"
plot(plotdata_zoom$Max_Blade_Length_m, plotdata_zoom$MaxStemDia_cm, col = col, cex=1, xlab="Maximum lamina length [m]", ylab="Maximum stem diameter [cm]")
abline(v=1.30)
abline(h=6)
dev.off()

#lot(traitData3$Max_Blade_Length_m, traitData3$Max_Petiole_length_m + traitData3$Max_Rachis_Length_m)

traitData3[!is.na(traitData3$MaxStemDia_cm) & traitData3$MaxStemDia_cm>6 & traitData3$NovGen == "Dypsis",c("figurenames", "MaxStemDia_cm", "Max_Blade_Length_m")]
traitData3[!is.na(traitData3$Max_Blade_Length_m) & traitData3$Max_Blade_Length_m>1.3 & traitData3$NovGen == "Dypsis",c("figurenames", "MaxStemDia_cm", "Max_Blade_Length_m")]

#plotdata[!is.na(plotdata$MaxStemDia_cm>5) & plotdata$MaxStemDia_cm>5 & plotdata$NovGen == "Dypsis",c( "MaxStemDia_cm", "Max_Blade_Length_m")]

traitData3[!is.na(traitData3$MaxStemDia_cm) & traitData3$MaxStemDia_cm<=6 & !is.na(traitData3$Max_Blade_Length_m) & traitData3$Max_Blade_Length_m <1.3 & traitData3$NovGen == "Chrysalidocarpus",c("figurenames", "MaxStemDia_cm", "Max_Blade_Length_m")]

# Alternative: leaf length
##########################

plotdata <- traitData3[traitData3$NovGen %in% c("Dypsis", "Chrysalidocarpus"),]
plotdata$NovGen <- as.factor(plotdata$NovGen)
#plotdata[plotdata$figurenames == "D. makirae", "Max_Blade_Length_m"] <- 0.15 + 0.25 # Petiole + Rachis, Rakotoarinivo et al. 2009 PALMS

plotdata[is.na(plotdata$Max_Blade_Length_m),]
plotdata <- plotdata[,c("NovGen", "MaxStemHeight_m", "MaxStemDia_cm", "MaxLeafNumber", "Max_Blade_Length_m", "Max_Petiole_length_m", "Max_Rachis_Length_m", "AverageFruitLength_cm", "AverageFruitWidth_cm"
                        , "MaxPeduncleLength_cm", "MaxRachisLength_cm")]

plotdata$MaxLeafLength <- plotdata$Max_Blade_Length_m + plotdata$Max_Petiole_length_m
plotdata$MaxPetioleRachis <- plotdata$Max_Petiole_length_m + plotdata$Max_Rachis_Length_m

plotdata <- plotdata[,colnames(plotdata) != c("Max_Rachis_Length_m")]
plotdata <- plotdata[,!colnames(plotdata) %in% c("Max_Blade_Length_m", "MaxPetioleRachis")]
#plotdata <- plotdata[,!colnames(plotdata) %in% c("Max_Blade_Length_m", "MaxLeafLength", "Max_Petiole_length_m")]

pdf(paste(figurepath, "/biplot_leaflength.pdf", sep=""), width=1.5*6.93, height = 1.5*6.93/2)
par(mfrow=c(1,2))

plotdata_zoom <- plotdata[plotdata$MaxStemDia_cm<10 & plotdata$MaxLeafLength<2.5,]

col = rep("black", nrow(plotdata))
col[plotdata$NovGen == "Dypsis"] <- "red"
plot(plotdata$MaxLeafLength, plotdata$MaxStemDia_cm, col = col, cex=1, xlab="Maximum leaf length [m]", ylab="Maximum stem diameter [cm]")
abline(v=1.65)
abline(h=6)
lines(c(-0.1,max(plotdata_zoom$MaxLeafLength, na.rm=T)), c(max(plotdata_zoom$MaxStemDia_cm, na.rm=T),max(plotdata_zoom$MaxStemDia_cm, na.rm=T)), lty=3)
lines(c(max(plotdata_zoom$MaxLeafLength, na.rm=T),max(plotdata_zoom$MaxLeafLength, na.rm=T)), c(max(plotdata_zoom$MaxStemDia_cm, na.rm=T),-2), lty=3)

col = rep("black", nrow(plotdata_zoom))
col[plotdata_zoom$NovGen == "Dypsis"] <- "red"
plot(plotdata_zoom$MaxLeafLength, plotdata_zoom$MaxStemDia_cm, col = col, cex=1, xlab="Maximum leaf length [m]", ylab="Maximum stem diameter [cm]")
abline(v=1.75)
abline(h=6)
dev.off()

traitData3[!is.na(traitData3$MaxStemDia_cm) & traitData3$MaxStemDia_cm>6 & traitData3$NovGen == "Dypsis",c("figurenames", "MaxStemDia_cm", "Max_Blade_Length_m", "Max_Petiole_length_m")]
traitData3[traitData3$Max_Petiole_length_m + traitData3$Max_Blade_Length_m > 1.75 & traitData3$NovGen == "Dypsis",c("figurenames", "MaxStemDia_cm", "Max_Blade_Length_m", "Max_Petiole_length_m")]

#plotdata[!is.na(plotdata$MaxStemDia_cm>5) & plotdata$MaxStemDia_cm>5 & plotdata$NovGen == "Dypsis",c( "MaxStemDia_cm", "MaxLeafLength")]

traitData3[!is.na(traitData3$MaxStemDia_cm) & traitData3$MaxStemDia_cm<=6 & !is.na(traitData3$MaxLeafLength) & traitData3$MaxLeafLength <1.3 & traitData3$NovGen == "Chrysalidocarpus",c("figurenames", "MaxStemDia_cm", "MaxLeafLength")]

outlier_dypsis <- c("D. coursii", "D. rivularis", "D. rosea", "D. pinnatifrons", "D. marojejyi", "D. brevicaulis")
traitData3[traitData3$NovGen == "Chrysalidocarpus" | traitData3$figurenames %in% outlier_dypsis,"Acaulescent"]
