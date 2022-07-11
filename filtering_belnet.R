library(data.table)
library(dplyr)
library(pbapply)
library(base)

#rules:
#-at least 10 reads per junction (for each gene, IJC + SJC per cell > 10)
#-AS event has to be detected in >10 cells (for each gene, IJC + SJC > 0 for at least 10 cells)
#-calculated PSI value = IJC/(IJC+SJC)
#-AS event cannot have no variability across cells (for each gene, variance of PSI values not equal to 0)

#input: table with rows as cells and IJC and SJC as columns
data <- fread("/Users/catherinez/Research/ISMB/belnet_IJC_SJC_396.csv")
data1 <- data[,-c(1:5)]

datam <- data.matrix(data1)
rownames(datam) <- data$barcode
colnames(datam) <- colnames(data1)

t.data <- t(datam)
write.csv(t.data, "/Users/catherinez/Research/ISMB/belnet_IJC_SJC_t.csv")

exonslist <- pblapply(seq(1,nrow(t.data)-1,by=2), function(i) {
  x <- rownames(t.data)
  final <- as.data.frame(rbind(t.data[i, ], t.data[(i+1), ]))
  rownames(final) <- c(x[i], x[i+1])
  return(final)
})

names(exonslist) <- unique(sub("_SJC|_IJC", "", rownames(t.data)))

#remove low junction reads
removelow <- pblapply(exonslist, function(df) {
  df[, colSums(df) > 10]
})

transpose <- pblapply(removelow, function(df) {
  x = t(as.matrix(df))
  as.data.frame(x)
})

filteredcols <- transpose[sapply(transpose, nrow)>0]
toremove <- transpose[sapply(transpose, nrow)<1]

names <- names(toremove)
exonslist.f1 <- exonslist
exonslist.f1[names] <- NULL

#remove low expression in cells
removelowexp <- pblapply(exonslist.f1, function(df) {
  df[, colSums(df) > 0]
})

transpose1 <- pblapply(removelowexp, function(df) {
  x = t(as.matrix(df))
  as.data.frame(x)
})

filteredcols1 <- transpose1[sapply(transpose1, nrow)>10]
toremove1 <- transpose1[sapply(transpose1, nrow)<11]

names1 <- names(toremove1)
exonslist.f2 <- exonslist.f1
exonslist.f2[names1] <- NULL

exonslist.f2.df <- do.call(rbind, exonslist.f2)
t.exonslist.f2.df <- t(exonslist.f2.df)
write.csv(exonslist.f2.df, "/Users/catherinez/Research/ISMB/belnet_IJC_SJC_filtered.csv")

#psi values
exons.psi <- pblapply(exonslist.f2, function(df) {
  df <- rbind ( df, df[1, ] / colSums(df) )
  df[-c(1,2),]
})

exons.psi.df <- do.call(rbind, exons.psi)
exons.psi.df[is.na(exons.psi.df)] = 0
write.csv(exons.psi.df, "/Users/catherinez/Research/ISMB/belnet_psi_values_filtered_raw.csv")

#remove 0 variance exons
exons.psi.filter_var <- exons.psi.df[apply(exons.psi.df, 1, var) != 0, ]
write.csv(exons.psi.filter_var, "/Users/catherinez/Research/ISMB/belnet_psi_values_filtered_var.csv")


exons.psi.filter_var.m <- data.matrix(exons.psi.filter_var)
rownames(datam) <- data$barcode
colnames(datam) <- colnames(data1)

t.exons.psi.filter_var.m <- t(exons.psi.filter_var.m)
write.csv(t.exons.psi.filter_var.m, "/Users/catherinez/Research/ISMB/belnet_psi_values_filtered_var_t.csv")

#create sums of IJC/SJC for subsets
meta <- data[,c(1:5)]
subset_list <- meta$subset
uniq_subs <- unique(subset_list[subset_list != ""])
names(uniq_subs) <- uniq_subs
datamdt <- data.table(datam)
sumsubs <- pblapply(uniq_subs, function(sub) {
  index <- grep(sub, subset_list)
  sum <- colSums(datamdt[index])
})
sumsubsm <- do.call(cbind, sumsubs)

write.csv(sumsubsm, "/Users/catherinez/Research/ISMB/belnet_subset_sums.csv")

#create sums of IJC/SJC for tumor/normal
tumor_status_list <- meta$type
tumor_status <- unique(tumor_status_list[tumor_status_list != ""])
names(tumor_status) <- tumor_status
sumtus <- pblapply(tumor_status, function(sub) {
  index <- grep(sub, tumor_status_list)
  sum <- colSums(datamdt[index])
})
sumtusm <- do.call(cbind, sumtus)

write.csv(sumtusm, "/Users/catherinez/Research/ISMB/belnet_tumor_status_sums.csv")

#test--------------------------------------------------------------------------
test.t.data <- t.data[c(29:50), ]

test.exonslist <- pblapply(seq(1,nrow(test.t.data)-1,by=2), function(i) {
  x <- rownames(test.t.data)
  final <- as.data.frame(rbind(test.t.data[i, ], test.t.data[(i+1), ]))
  rownames(final) <- c(x[i], x[i+1])
  return(final)
})

names(test.exonslist) <- unique(sub("_SJC|_IJC", "", rownames(test.t.data)))

#remove low junction reads
test.removelow <- pblapply(test.exonslist, function(df) {
  df[, colSums(df) > 10]
})

test.transpose <- pblapply(test.removelow, function(df) {
  x = t(as.matrix(df))
  as.data.frame(x)
})

test.filteredcols <- test.transpose[sapply(test.transpose, nrow)>0]
test.toremove <- test.transpose[sapply(test.transpose, nrow)<1]

names <- names(test.toremove)

test.exonslist.f1 <- test.exonslist
test.exonslist.f1[names] <- NULL

#remove low expression in cells

test.removelowexp <- pblapply(test.exonslist.f1, function(df) {
  df[, colSums(df) > 0]
})

test.transpose1 <- pblapply(test.removelowexp, function(df) {
  x = t(as.matrix(df))
  as.data.frame(x)
})

test.filteredcols1 <- test.transpose[sapply(test.transpose, nrow)>10]

test.toremove1 <- test.transpose1[sapply(test.transpose1, nrow)<11]

test.names1 <- names(test.toremove1)

test.exonslist.f2 <- test.exonslist.f1
test.exonslist.f2[test.names1] <- NULL

test.exonslist.f2

test.exons.psi <- pblapply(test.exonslist.f2, function(df) {
  df <- rbind ( df, df[1, ] / colSums(df) )
  df[-c(1,2),]
})

test.exons.psi.df <- do.call(rbind, test.exons.psi)
test.exons.psi.df[is.na(test.exons.psi.df)] = 0

test.psi.filter_var <- test.exons.psi.df[apply(test.exons.psi.df, 1, var) >= 0.0001, ]

  
