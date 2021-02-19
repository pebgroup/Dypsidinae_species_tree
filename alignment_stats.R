base_path <- "/Users/au265104/Documents/WOLF/PROJECTS/65 Dypsis systematics paper/~manuscript/stats"

# Set working directory
# For outgroup/subtribes analysis:
# wd <- "~subtribes/AMAS" 
# For main analysis:

# Set number of taxa for occupancy calculation
# Outgroup/subtribes
# N <- 19
# Main analysis:
N <- 173

wd <- "~main/AMAS"
setwd(paste(base_path,"/",wd, sep=""))

summary_all <- read.table("summary_all.txt", header=T)
summary_exon <- read.table("summary_exon.txt", header=T)
summary_intron <- read.table("summary_intron.txt", header=T)

stats <- matrix(nrow=20,ncol=6)
colnames(stats) <- c("VAR","MEDIAN","MEAN","SD","MIN","MAX")
stats <- as.data.frame(stats)
stats$VAR <- c("No_of_taxa","Occupancy",rep(c("Alignment_length","Missing_percent","No_variable_sites","Proportion_variable_sites", "Parsimony_informative_sites", "Proportion_parsimony_informative"),3))

# supercontigs

summary_all <- summary_all[,c("No_of_taxa","Alignment_length","Missing_percent","No_variable_sites","Proportion_variable_sites", "Parsimony_informative_sites", "Proportion_parsimony_informative")]
summary_all$Occupancy <- 100*summary_all$No_of_taxa/N
summary_all <- summary_all[,c("No_of_taxa","Occupancy","Alignment_length","Missing_percent","No_variable_sites","Proportion_variable_sites", "Parsimony_informative_sites", "Proportion_parsimony_informative")]
summary_all$Proportion_parsimony_informative <- 100*summary_all$Proportion_parsimony_informative
summary_all$Proportion_variable_sites <- 100*summary_all$Proportion_variable_sites

stats$MEDIAN[1:8] <-  sapply(summary_all,median)
stats$MEAN[1:8] <-  sapply(summary_all,mean)
stats$SD[1:8] <-  sapply(summary_all,sd)
stats$MIN[1:8] <-  sapply(summary_all,min)
stats$MAX[1:8] <-  sapply(summary_all,max)

# exons

summary_exon <- summary_exon[,c("Alignment_length","Missing_percent","No_variable_sites","Proportion_variable_sites", "Parsimony_informative_sites", "Proportion_parsimony_informative")]

summary_exon$Proportion_parsimony_informative <- 100*summary_exon$Proportion_parsimony_informative
summary_exon$Proportion_variable_sites <- 100*summary_exon$Proportion_variable_sites

stats$MEDIAN[9:14] <-  sapply(summary_exon,median)
stats$MEAN[9:14] <-  sapply(summary_exon,mean)
stats$SD[9:14] <-  sapply(summary_exon,sd)
stats$MIN[9:14] <-  sapply(summary_exon,min)
stats$MAX[9:14] <-  sapply(summary_exon,max)

# introns

summary_intron <- summary_intron[,c("Alignment_length","Missing_percent","No_variable_sites","Proportion_variable_sites", "Parsimony_informative_sites", "Proportion_parsimony_informative")]

summary_intron$Proportion_parsimony_informative <- 100*summary_intron$Proportion_parsimony_informative
summary_intron$Proportion_variable_sites <- 100*summary_intron$Proportion_variable_sites

stats$MEDIAN[15:20] <-  sapply(summary_intron,median)
stats$MEAN[15:20] <-  sapply(summary_intron,mean)
stats$SD[15:20] <-  sapply(summary_intron,sd)
stats$MIN[15:20] <-  sapply(summary_intron,min)
stats$MAX[15:20] <-  sapply(summary_intron,max)

stats[,2:6] <- round(stats[,2:6],2)

write.table(stats, "alignment_stats.csv", dec = ",", sep=";")

sum(summary_all$No_of_taxa)



