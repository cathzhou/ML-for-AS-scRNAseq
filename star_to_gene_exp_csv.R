library(reshape2)
library(tibble)
library(gtools)
library(dplyr)

setwd("/oak/stanford/groups/ebutcher/catherine/PLN_EC/star/output")

format <- function(id) {
  cell <- read.delim(paste0(getwd(),"/",id,"ReadsPerGene.out.tab"))
  cell <- cell[-c(1:3),]
  names(cell)[2] <- id
  cell <- cell[ -c(1,3,4) ]
}

cellids = readLines("/oak/stanford/groups/ebutcher/catherine/PLN_EC/csv/cellids.txt")
allcells = lapply(cellids, format)
master <- read.delim("/oak/stanford/groups/ebutcher/catherine/PLN_EC/csv/belnet_genelist.txt")
allcells <- append(allcells, value="0", after = 0)
allcells[[1]] <- data.frame(master)
final = do.call(cbind, allcells)

write.csv(final, "/oak/stanford/groups/ebutcher/catherine/PLN_EC/csv/raw_gene_expression.csv")

rm(list=ls())