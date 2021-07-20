library(phytools)
library(writexl)
library(geiger)
library(dplyr)
library(stringr)

setwd("~/OneDrive - Aarhus Universitet/PROJECTS/65 Dypsis systematics paper/~manuscript/figures/trait figure")

traitData3 <- read.csv2("DypsidinaeTraitData.csv",header = TRUE)
traitData3 <- traitData3[!traitData3$SpecName=="",]
traitData3$SpecName[traitData3$SpecName == "Dypsis Leucomalla"] <- "Dypsis leucomalla"
traitData3 <- traitData3[!traitData3$SpecName=="Loxococcus rupicola",] #exclude Loxococcus

read.tree("astral_tree.tre") -> tree
tree$edge.length[is.na(tree$edge.length)] <- 0.001
write.tree(tree,"astral_tree_to_ladderize.tre") #ladderize in FigTree - ladderize() messes up the plot
read.tree("astral_tree_lad.tre") -> tree
tree <- drop.tip(tree, c("0194", "0196", "0199", "0202", "0204", "1012", "1011"))

# match tree tips to trait dataset
##################################

read.table("rename.csv", sep=";", colClasses = "character") -> rename
rownames(rename) <- rename$V1

names <- rename[tree$tip.label,"V2"]
names[names == "Dypsis-heteromorpha-JD7822_S62_L001"] = "Dypsis-sp-heteromorpha-JD7822_S62_L001"
names[names == "Dypsis-mirabilis-florencei_S64_L001"] = "Dypsis-sp-mirabilis-florencei_S64_L001"
names[names == "Dypsis-onilahensis-Isalo-form_S56_L001"] = "Dypsis-sp-onilahensis-Isalo-form_S56_L001"
names <- str_replace_all(names,"-","_")
names2 <- c()
for(i in strsplit(names,"_")){
  names2 <- c(names2, paste(i[1:2], collapse="_"))
}
names2[names2 == "Dypsis_sp"] <- paste("Dypsis_sp", 1:8, sep="_")

read.table("figurenames2.csv", sep=";", colClasses = "character")[1:174,] -> rename2 
rownames(rename2) <- rename2$V1
figurenames <- rename2[tree$tip.label,"V2"]
names(figurenames) <- names2

traitData3$figurenames <- figurenames[str_replace(traitData3$SpecName," ","_")]

# rename tree with figure names
###############################
tree$tip.label <- figurenames

traitData3$figurenames

# Build plot
############

pdf("traitplot.pdf", width=8.3, height = 11.7)

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

legend("topleft", c("Max. stem height", "Max. stem diameter", "Max. blade length", "Fruit volume", "Peduncle length", "Inflorescence length", "3 stamens", "6 stamens", "52-59 stamens", "Antesepalous stamens", "Antepetalous stamens", "Smooth endosperm", "Ruminate endosperm"), 
       col = c("brown", "blue","purple","blue4", "forestgreen", "red", "orange", "orange", "orange", "turquoise3", "turquoise3", "grey31", "grey31"),
       text.col = "black", pch = c(19,19,19,19,19,19,1,19,8,1,19,1,19), bg = "white")

dev.off()


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


#traitData3$SpecName <- str_replace(traitData3$SpecName," ","_")
#name.check(tree,traitData3$SpecName)
############ continuous traits #########
#########Max stem hieght####
traitData3 <- traitData3[1:155,]
rownames(traitData3) <- traitData3$SpecName[1:155]
stem.height <- traitData3$MaxStemHeight_m
log.stem <- log10(traitData3$MaxStemHeight_m)
hist(log.stem)
names(log.stem) <- str_replace(traitData3$SpecName," ","_")
foo<-function(tree,data){
  fsize<-36*par()$pin[2]/par()$pin[1]/Ntip(tree)
  plotTree.barplot(tree,data,args.plotTree=list(fsize=fsize),lwd=1)
}
pdf("stemheight.pdf", height=11.75, width=8.25)
foo(tree,log.stem)
dev.off()

####MAx Stem Diameter###
dia<- traitData3$MaxStemDia_cm
names(dia) <- str_replace(traitData3$SpecName," ","_")
pdf("dia.barplot.pdf", height=11.75, width=8.25)
foo(tree,dia)
dev.off()

pdf("Dypsis_tree2.pdf", height=11.75, width=8.25)
plotTree(new.tree)
dev.off()

#####Max Blade length###
blade_length<- traitData3$Max_Blade_Length_m
names(blade_length) <- str_replace(traitData3$SpecName," ","_")
pdf("blade_length.barplot.pdf", height=11.75, width=8.25)
foo(tree,blade_length)
dev.off()

####Max rachis length###
rachis <- traitData3$Max_Rachis_Length_m
names(rachis) <- str_replace(traitData3$SpecName," ","_")
pdf("Rachis.pdf", height=11.75, width=8.25)
foo(tree,rachis)
dev.off()

#####Max petiole lengh######
petiole <- traitData3$Max_Petiole_length_m
names(petiole) <- str_replace(traitData3$SpecName," ","_")
pdf("MAx_Petiole_length.pdf", height=11.75, width=8.25)
foo(tree,petiole)
dev.off()

#####Max petiole lengh######
petiole <- traitData3$Max_Petiole_length_m
names(petiole) <- str_replace(traitData3$SpecName," ","_")
pdf("MAx_Petiole_length.pdf", height=11.75, width=8.25)
foo(tree,petiole)
dev.off()

####MAx leaf number###
max_leaf_n <- traitData3$MaxLeafNumber
names(max_leaf_n) <- str_replace(traitData3$SpecName," ","_")
pdf("max_leaf_n.barplot.pdf", height=11.75, width=8.25)
foo(tree,max_leaf_n)
dev.off()

####Average fruit length###
traitData3$AverageFruitLength_cm
fruit_length<- traitData3$AverageFruitLength_cm
names(fruit_length) <- str_replace(traitData3$SpecName," ","_")
pdf("Fruit_length.barplot.pdf", height=11.75, width=8.25)
foo(tree,fruit_length)
dev.off()

####Average fruit width####
fruit_width <- traitData3$AverageFruitWidth_cm
names(fruit_width) <- str_replace(traitData3$SpecName," ","_")
pdf("Fruit_width.barplot.pdf", height=11.75, width=8.25)
foo(tree,fruit_width)
dev.off()

###Peduncle-Inflorescence Length Proportion:
PedunclePro <- traitData3$PeduncleInflorescenceLengthProportion
names(PedunclePro) <- str_replace(traitData3$SpecName," ","_")
pdf("PedunInfloProper.barplot.pdf", height=11.75, width=8.25)
foo(tree,PedunclePro)
dev.off()


#####Let's try categorical variables####
canopy <- traitData3$UnderstoreyCanopy
names(canopy) <- str_replace(traitData3$SpecName," ","_")
pdf("Understory_Canopy.pdf", height=11.75, width=8.25)
foo1<-function(tree,data){
  fsize<-36*par()$pin[2]/par()$pin[1]/Ntip(tree)
  dotTree(tree,data,fsize=0.4,colors=setNames(c("blue","red"),c("understorey","canopy")),ftype="i",lwd=1)
}

foo1(tree,canopy)
dev.off()

##### Ruminate /homogeneous endosperm:
endosperm <- as.character(traitData3$RuminateEndosperm)
names(endosperm) <- str_replace(traitData3$SpecName," ","_")
pdf("Endosperm.pdf", height=11.75, width=8.25)
foo1<-function(tree,data){
  fsize<-36*par()$pin[2]/par()$pin[1]/Ntip(tree)
  dotTree(tree,data,fsize=0.4,colors=setNames(c("blue","red"),c("0","1")),ftype="i",lwd=1)
}

foo1(tree,endosperm)
dev.off()

dotTree(tree,endosperm,fsize=0.2,colors=setNames(c("blue","red"),c("0","1")),ftype="i",lwd=1,border=0.1, fsize=0.1,outline=TRUE)
dotTree(tree,canopy,fsize=0.2,colors=setNames(c("blue","red"),c("0","1")),ftype="i",lwd=1,border=0.1, fsize=0.1,outline=TRUE,add=TRUE)

#### Stamen number:
stamen <- as.character(traitData3$StamenNumber)
names(stamen) <- str_replace(traitData3$SpecName," ","_")
pdf("Stamen_Number.pdf", height=11.75, width=8.25)
foo1<-function(tree,data){
  fsize<-36*par()$pin[2]/par()$pin[1]/Ntip(tree)
  dotTree(tree,data,fsize=0.4,colors=setNames(c("blue","red"),c("3","6")),ftype="i",lwd=1)
}

foo1(tree,stamen)
dev.off()

### Antepetalous
Antepetalous <- as.character(traitData3$Antepetalous)
names(Antepetalous) <- str_replace(traitData3$SpecName," ","_")
pdf("Antepetalous.pdf", height=11.75, width=8.25)
foo1<-function(tree,data){
  fsize<-36*par()$pin[2]/par()$pin[1]/Ntip(tree)
  dotTree(tree,data,fsize=0.4,colors=setNames(c("blue","red"),c("0","1")),ftype="i",lwd=1)
}

foo1(tree,Antepetalous)
dev.off()

###Crownshaft presence
crownshaft <- as.character(traitData3$Crownshaft)
names(crownshaft) <- str_replace(traitData3$SpecName," ","_")
pdf("Crownshaft.pdf", height=11.75, width=8.25)
foo1<-function(tree,data){
  fsize<-36*par()$pin[2]/par()$pin[1]/Ntip(tree)
  dotTree(tree,data,fsize=0.4,colors=setNames(c("blue","red"),c("0","1")),ftype="i",lwd=1)
}

foo1(tree,crownshaft)
dev.off()

#Anther form:
antherform <- as.character(traitData3$AntherForm)
names(antherform) <- str_replace(traitData3$SpecName," ","_")
pdf("Antherform.pdf", height=11.75, width=8.25)
foo1<-function(tree,data){
  fsize<-36*par()$pin[2]/par()$pin[1]/Ntip(tree)
  dotTree(tree,data,fsize=0.4,colors=setNames(c("blue","red","yellow","green","purple","brown","orange"),c("didymous","obtuse","elongate","sagittate","pedunculous","cornical","globose")),ftype="i",lwd=1)
}

foo1(tree,antherform)
dev.off()

######## Inflorescence branching order
branhingorder <- as.character(traitData3$InflorescenceBranchingOrder)
names(branhingorder) <- str_replace(traitData3$SpecName," ","_")
pdf("InfloBranchingOrder.pdf", height=11.75, width=8.25)
foo1<-function(tree,data){
  fsize<-36*par()$pin[2]/par()$pin[1]/Ntip(tree)
  dotTree(tree,data,fsize=0.4,colors=setNames(c("blue","red","yellow","green","purple","brown","orange","pink"),c("0","1","1-2","2","2-3","3","3-4","4")),ftype="i",lwd=1)
}

foo1(tree,branhingorder)
dev.off()

#### Stem solitary or not###
solitary <- traitData3$StemSolitary
names(solitary) <- str_replace(traitData3$SpecName," ","_")
pdf("Solitary_Canopy.pdf", height=11.75, width=8.25)
foo1<-function(tree,data){
  fsize<-36*par()$pin[2]/par()$pin[1]/Ntip(tree)
  dotTree(tree,data,fsize=0.5,colors=setNames(c("blue","red","green"),c("0","1","2")),ftype="i",lwd=1)
}

foo1(tree,solitary)
dev.off()

#### fruit shape #######
fruit_shape <- traitData3$FruitShape
names(fruit_shape) <- str_replace(traitData3$SpecName," ","_")
pdf("fruit_shape.pdf", height=11.75, width=8.25)
foo1<-function(tree,data){
  fsize<-36*par()$pin[2]/par()$pin[1]/Ntip(tree)
  dotTree(tree,data,fsize=0.5,colors=setNames(c("blue","red","green","yellow"),c("ovoid","globose","elongate","ellipsoid")),ftype="i",lwd=1)
}

foo1(tree,fruit_shape)
dev.off()

######Fruit colour#####
#20 different fruit colours - I simplify the dataset and only takes the first mentioned colour for each species###
col_names<- strsplit(traitData3$MainFruitColors, "\\;")
h1 <- as.list("")
for(i in 1:173) {
  h1[[i]] <- col_names[[i]] [1]
  print(h1[i])
}
h1 <- unlist(h1)
h1 <- as.character(h1)
names(h1) <- str_replace(traitData3$SpecName," ","_")
pdf("fruit_col.pdf", height=11.75, width=8.25)
foo1<-function(tree,data){
  fsize<-36*par()$pin[2]/par()$pin[1]/Ntip(tree)
  dotTree(tree,data,fsize=0.5,colors=setNames(c("red","orange","green","yellow","black","purple","brown"),c("red","orange","green","yellow","black","purple","brown")),ftype="i",lwd=1)
}

foo1(tree,h1)
dev.off()

################Phylomorphospace of PCs 1 and 2:
#install.packages("dplyr")
library(dplyr)
traitData3$SpecName <- str_replace(traitData3$SpecName," ","_")
library(dplyr)
sizes <- traitData3 %>% select(MaxStemHeight_m, MaxStemDia_cm, Max_Blade_Length_m)
rownames(sizes) <- traitData3$SpecName
sizes <- na.omit(sizes)
sizes$MaxStemHeight_m
library(geiger)

chk<-name.check(ecomorph.tree,sizes)
chk

#Dypsis_sanctaemariae is placed twice in phylogeny.. Only keep one:
unique_names<-unique(ecomorph.tree$tip.label)
ecomorph.tree <- keep.tip(ecomorph.tree, unique_names)

#Only include species that is both in dataset and tree:
ecomorph.tree<-drop.tip(tree,chk$tree_not_data)
name.check(ecomorph.tree,sizes)
#difine the six clades:
plotTree(ecomorph.tree,fsize=0.3,ftype="i",lwd=0.5)
nodelabels(bg="white",cex=0.3,frame = "circle")
tiplabels(bg="white",cex=0.3,frame = "circle")

chrysalidocarpus<-extract.clade(ecomorph.tree,209)
dypsis<-extract.clade(ecomorph.tree,139)
marojejya<-extract.clade(ecomorph.tree,137)
vonitra<-extract.clade(ecomorph.tree,129)
lemurophoenix<-as.character("Lemurophoenix_halleuxii")
masoala<-extract.clade(ecomorph.tree,126)
loxococcus<- "Loxococcus_rupicola"
species.list <- (c(loxococcus, masoala$tip.label,lemurophoenix,vonitra$tip.label,
                   marojejya$tip.label,dypsis$tip.label,chrysalidocarpus$tip.label))
phylo.group<- as.vector(length(ecomorph.tree$tip.label))
phylo.group[1]<-"Loxococcus"
phylo.group[2:3]<-"Masoala"
phylo.group[4]<-"Lemurophoenix"
phylo.group[5:12]<-"Vonitra"
phylo.group[13:14]<-"Marojejya"
phylo.group[15:85]<-"Dypsis"
phylo.group[86:123]<-"Chrysalidocarpus"

groups <- as.data.frame(phylo.group)
rownames(groups)<-species.list

#The data should be log transformed as the data in sizes have different orders of magnitude:
sizes$MaxStemHeight_m <- log10(sizes$MaxStemHeight_m)
sizes$MaxStemDia_cm <- log10(sizes$MaxStemDia_cm)
sizes$Max_Blade_Length_m <- log10(sizes$Max_Blade_Length_m)
sizes <- na.omit(sizes)


#Run phylo PCA 
ecomorph.pca<-phyl.pca(ecomorph.tree,sizes)
str(ecomorph.pca)
scores(ecomorph.pca)[,1:3]

#See how much the different PC axis explains:
par(mar=c(4.1,4.1,2.1,1.1),las=1)
plot(ecomorph.pca,main="")

#Plot the first two PC axis:
plot(scores(ecomorph.pca)[,1:2])

#Make PhyloPCA plot:
phylomorphospace(ecomorph.tree,scores(ecomorph.pca)[,1:2],node.size=c(0,0.3),bty="n",ftype="off")
eco<-setNames(groups[,1],rownames(groups))
eco<-as.factor(eco)
ECO<-to.matrix(eco,levels(eco))
tiplabels(pie=ECO[ecomorph.tree$tip.label,],cex=0.4)
legend(x="bottomleft",legend=levels(eco),cex=0.8,pch=21,
       pt.bg=rainbow(n=length(levels(eco))),pt.cex=1.5)

#plot phyloPCA with species names:
eco<-setNames(nm = rownames(sizes))

eco<-as.factor(eco)

project.phylomorphospace(ecomorph.tree,fsize=0.4, scores(ecomorph.pca)[,1:2],node.size=c(0,0.05),pie=ECO[ecomorph.tree$tip.label,],cex=0.1)

ECO<-to.matrix(eco,levels(eco))

pdf("Phylo_PCA.pdf", height=10, width=10)
plotTree(tree1)
dev.off()

#Barplots with the scores from the two first axis in the phylo PCA:
PC1.score <- setNames(scores(ecomorph.pca)[,1],rownames(ecomorph.pca$S))
foo2<-function(tree,data){
  fsize<-36*par()$pin[2]/par()$pin[1]/Ntip(tree)
  plotTree.barplot(tree,data,args.plotTree=list(fsize=fsize),lwd=1,xlim=c(-80,20))
}
pdf("Barplot_PC1.scores.pdf", height=11.75, width=8.25)
foo(ecomorph.tree,PC1.score)
dev.off()

#######PC2 #######
PC2.score <- setNames(scores(ecomorph.pca)[,2],rownames(ecomorph.pca$S))
foo2<-function(tree,data){
  fsize<-36*par()$pin[2]/par()$pin[1]/Ntip(tree)
  plotTree.barplot(tree,data,args.plotTree=list(fsize=fsize),lwd=1)
}
pdf("Barplot_PC2.scores.pdf", height=11.75, width=8.25)
foo(ecomorph.tree,PC2.score)
dev.off()

