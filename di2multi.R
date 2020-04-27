library(ape)

args <- commandArgs(trailingOnly=TRUE)

path <- args[1]

tree <- read.tree(path)

tree <- di2multi(tree, tol=0.00001)

write.tree(tree, paste(path, ".noshort", sep=""))

print(path)