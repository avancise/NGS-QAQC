rm(list = ls())
source("callBasesPileup1.R")

setwd("C:/Users/amv/Google Drive/11 Global Taxonomy/Mitogenomes/QAQC")
plp.fname <- "Gmac_mito_mtNDA_lib2.pileup.ref.csv"
x <- callBasesPileup(plp.fname)


seq.mat <- as.character(x$dna.seqs)
#remove failed sequences and padded ends
seq.mat<-seq.mat[!rownames(seq.mat) %in% c("z0157267","z0175342","z0175343","z0175843","z0175957"),]
seq.mat<-seq.mat[,-seq(from=1,to=40,1)]
seq.mat<-seq.mat[,-seq(from=16388,to=16427,1)]

#count number of Ns
num.ns <- apply(seq.mat, 1, function(z) sum(tolower(z) == "n"))
sum(num.ns)

#remove sequences with >10 Ns
high.ns<-num.ns[which(num.ns>10)]
high.n.names<-names(high.ns)
seq.mat<-seq.mat[!rownames(seq.mat) %in% high.n.names,]

#recount number of Ns
num.ns.2<-apply(seq.mat, 1, function(z) sum(tolower(z) == "n"))
sum(num.ns.2)
ns<-num.ns.2[which(num.ns>0)]

#save edited file as .fasta
fname <- gsub(".csv", "edit4.fasta", plp.fname)
write.dna(
  seq.mat, fname, format = "fasta", nbcol = -1,
  colsep = "", indent = 0, blocksep = 0
)

#matrix to dataframe
library(reshape2)
seq.df<-setNames(melt(seq.mat), c('ID','position','allele'))
seq.df$allele<-as.character(seq.df$allele)
seq.df.ns<-seq.df[which(seq.df$allele=="n"),]
