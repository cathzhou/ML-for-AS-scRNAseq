.libPaths("~/R/x86_64-pc-linux-gnu-library/4.0")
library(reshape2)
library(tibble)
library(gtools)
library(dplyr)

setwd("/oak/stanford/groups/ebutcher/catherine/PLN_EC/rmats/output/merged_output")

format <- function(id) {
  cell <- read.delim(paste0(getwd(),"/",id,"/SE.MATS.JC.txt"))
  cell <- cell[ -c(1,2,8:12,15:23) ]
  
  names(cell)[6] <- 'IJC'
  names(cell)[7] <- 'SJC'
  
  cell <- melt(cell, id.vars=1:5, value.name = id)
  cell$ident <- "0"
  cell$ident <- paste(cell$geneSymbol, cell$variable, 
                              cell$chr, cell$strand, cell$exonStart_0base,
                              cell$exonEnd, sep = "_")
  cell <- cell[-c(1:6)]
  cell <- cell[c(2,1)]
  cell %>% distinct(ident, .keep_all = TRUE)
}

loj <- function(x, y) {
  df <- merge(x = x, y = y, by = "ident", all.y = TRUE)
  df <- df[-c(1)]
}

cellids = readLines("/oak/stanford/groups/ebutcher/catherine/PLN_EC/csv/cellids.txt")
allcells = lapply(cellids, format)
master <- read.delim("/oak/stanford/groups/ebutcher/catherine/PLN_EC/csv/merged_csv/finaluniqexons.txt")
allcells = lapply(allcells, loj, master)
allcells <- append(allcells, value="0", after = 0)
allcells[[1]] <- data.frame(master)
final = do.call(cbind, allcells)

write.csv(final, "/oak/stanford/groups/ebutcher/catherine/PLN_EC/csv/merged_csv/finalallcells.csv")

rm(list=ls())

