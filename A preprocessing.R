##opening libraries, defining functions, and preprocessing

library(topGO)
library(ggplot2)
library(dplyr)
library(gtools)
library(moments)
library(tidyverse)
library(GenomicFeatures)
library(ggrepel)
library(ggplot2)
library(ggpmisc)
library(enrichR)

get_all_symbols = function(ID) {
  symbols_list<-unique(allcompgb2$Gene_Symbol[which(allcompgb2$TransID %in% ID)])
  toString(as.list(symbols_list))
}

symbols2ids <- function(symbols) {
  print("Warning: IDs are not in same order!")
  symbolIDs$id[which(symbolIDs$symbol %in% symbols)]
}

nthroot = function(x,n) {
  (abs(x)^(1/n))*sign(x)
}

setwd("C:/Users/jardi/OneDrive/Documents/TapinosLab/singlenucm6a")
##read in 
allcompgb2<-read.csv("C:/Users/jardi/OneDrive/Documents/TapinosLab/singlenucm6a/Report/PaperData/diffmethgenes/allcomparison/GB2ED.csv")
##fix colname formatting
colnames(allcompgb2)[which(names(allcompgb2)=="ï..TransID")]<-"TransID"
colnames(allcompgb2)[11:12] <- c("E.Mean.pm6A", "D.Mean.pm6A")
colnames(allcompgb2)[17:20] <- c("E1.pm6A", "E2.pm6A", "D1.pm6A", "D2.pm6A")
##convert percentage characters to decimal
##skip these and do them in Excel, otherwise Excel inserts weird characters into the csv
# allcompgb2$D1.pm6A<-as.numeric(sub("%", "", allcompgb2$D1.pm6A))
# allcompgb2$D2.pm6A<-as.numeric(sub("%", "", allcompgb2$D2.pm6A))
# allcompgb2$E1.pm6A<-as.numeric(sub("%", "", allcompgb2$E1.pm6A))
# allcompgb2$E2.pm6A<-as.numeric(sub("%", "", allcompgb2$E2.pm6A))
# allcompgb2$E.Mean.pm6A<-as.numeric(sub("%", "", allcompgb2$E.Mean.pm6A))
# allcompgb2$D.Mean.pm6A<-as.numeric(sub("%", "", allcompgb2$D.Mean.pm6A))

#normalize: divide m6A percent abundance by largest value,
#because having 2400% abundance makes no sense esp if we want to track entropy
maxenrichment<-max(allcompgb2[17:20]) #for global normalization
#sum(allcompgb2$E1.pm6A>1) #1504
#sum(allcompgb2$E2.pm6A>1) #567
#sum(allcompgb2$D1.pm6A>1) #1228
#sum(allcompgb2$D2.pm6A>1) #591
#maxd<-max(max(allcompgb2$D1.pm6A),max(allcompgb2$D2.pm6A))
#maxe<-max(max(allcompgb2$E1.pm6A),max(allcompgb2$E2.pm6A))
#max1<-max(max(allcompgb2$E1.pm6A),max(allcompgb2$D1.pm6A))
#max2<-max(max(allcompgb2$E2.pm6A),max(allcompgb2$D2.pm6A))

allcompgb2$D1.pm6A<-allcompgb2$D1.pm6A/maxenrichment#max1#100#max(allcompgb2[17:20])##max(allcompgb2$D1.pm6A) #these for batch normalization
#max2<-max(allcompgb2[17:20])
allcompgb2$D2.pm6A<-allcompgb2$D2.pm6A/maxenrichment#max2#100#max(allcompgb2[17:20])#maxenrichment#max(allcompgb2$D2.pm6A)
#max3<-max(allcompgb2[17:20])
allcompgb2$E1.pm6A<-allcompgb2$E1.pm6A/maxenrichment#max1#100#max(allcompgb2[17:20])#maxenrichment#max(allcompgb2$E1.pm6A)
#max4<-max(allcompgb2[17:20]) #previously, this value was much lower than the others, since the absolute max occurs in E1
#thus in the original normalization procedure the entacapone entropies were artificially inflated
allcompgb2$E2.pm6A<-allcompgb2$E2.pm6A/maxenrichment#max2#max(allcompgb2[17:20])#maxenrichment#max(allcompgb2$E2.pm6A)
#normalization method for these doesn't matter as much since they are not used in calculations
allcompgb2$E.Mean.pm6A<-allcompgb2$E.Mean.pm6A/maxenrichment#max(allcompgb2$E.Mean.pm6A)
allcompgb2$D.Mean.pm6A<-allcompgb2$D.Mean.pm6A/maxenrichment#max(allcompgb2$D.Mean.pm6A)
allcompgb2$Foldchange<-log(allcompgb2$E.Mean.pm6A/allcompgb2$D.Mean.pm6A)

write.csv(allcompgb2, "allcompgb2.csv")