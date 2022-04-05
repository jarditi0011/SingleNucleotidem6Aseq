##Entropy calculation

IDs<-allcompgb2$TransID
siteslist<-table(factor(IDs, levels=unique(IDs)))
#siteslist<-unique(allcompgb2$Gene_Symbol)
siteslist<-data.frame(siteslist)
num.tx=nrow(siteslist)
colnames(siteslist)[1]<-"TransID"

#create new cols for siteslist
siteslist$TransID<-as.character(siteslist$TransID)
siteslist$E1Entropy<-0 #3
siteslist$E2Entropy<-0
siteslist$D1Entropy<-0
siteslist$D2Entropy<-0
#entropy / number of sites, allows comparison of txs
siteslist$Norm.E1Entropy<-0 #7
siteslist$Norm.E2Entropy<-0
siteslist$Norm.D1Entropy<-0
siteslist$Norm.D2Entropy<-0
#compiling that data
siteslist$E_AvgEnt<-0 #11
siteslist$D_AvgEnt<-0
siteslist$EntDiff<-0
#string list of ints
#indicating sites of modification for each tx,
#in ascending order
siteslist$SiteLoci<-"NULL" #TODO 14
#string list of states formatted as 001101, 
#where each digit is a site and 0 indicates unmodified
siteslist$States<-"NULL" #TODO 15
#string list of probabilities, 
#indices of which correspond to each state in States
siteslist$StateProbs<-"NULL" #TODO 16
#L2 Measure (Euclidean Distance)
siteslist$Euclidean<-0 #17
siteslist$Gene_Symbol<-sapply(siteslist$TransID, get_all_symbols)

starttime<-Sys.time()
print(paste("Started at ", starttime))
#determine the entropy for each transcript, as well as the most likely state
for (txnum in 1:num.tx) {
  #get all sites on this tx
  my.tx<-siteslist$TransID[txnum]
  print(paste("transcript ", as.character(txnum), " of ", as.character(num.tx)))
  my.tx.sites<-allcompgb2[which(allcompgb2$TransID==my.tx),]
  
  #order the sites by 5' --> 3'
  my.tx.sites<-my.tx.sites[order(my.tx.sites$m6A_transcript_location),]
  num.tx.sites<-nrow(my.tx.sites)
  siteslist$SiteLoci[txnum]<-toString(as.list(my.tx.sites$m6A_transcript_location))
  print(paste(as.character(my.tx), " has ", as.character(num.tx.sites), " sites"))
  
  #isolate the probabilities of each site from my.tx.sites
  probs.sites<-my.tx.sites[17:20]
  
  #siteslist$SiteProbs[txnum]<-toString(as.list(probs.sites))
  
  #calculate entropy for each transcript by summing site entropies
  #H(x1...xn)=H(x1)+...+H(xn)
  pe<-my.tx.sites[11]
  pd<-my.tx.sites[12]
  siteslist$SiteProbs[txnum]<-toString(as.list(my.tx.sites$Foldchange))
  siteslist[txnum,][17]<-sqrt(sum((pe-pd)^2)/num.tx.sites)*sign(sum(pe-pd))#sum(pe-pd)/num.tx.sites#sqrt(sum((pe-pd)^2)/num.tx.sites^2)
  for (site in num.tx.sites) {
    p<-probs.sites[site,]
    siteentropy<-(-1)*(p*log(p)+(1-p)*log(1-p))
    #siteDiffSquared<-(pe[txnum,1]-pd[txnum,1])^2
    siteslist[txnum,][3:6]<-siteslist[txnum,][3:6]+siteentropy
  }
  siteslist[is.na(siteslist)]<-0
  #normalize entropies to same number of sites, 
  #allowing comparisons between txs with diff numbers of sites
  siteslist[txnum,][7:10]<-siteslist[txnum,][3:6]/num.tx.sites
  #siteslist[txnum,][17]<-sqrt(siteslist[txnum,][17]/num.tx.sites)
}
endtime=Sys.time()
print(paste("finished at ", endtime))
print(paste("Time used to calculate entropy:", 
            difftime(endtime, starttime, 
                     units = "secs"), "secs"))

siteslist$E_AvgEnt<-(siteslist$Norm.E1Entropy+siteslist$Norm.E2Entropy)/2
siteslist$D_AvgEnt<-(siteslist$Norm.D1Entropy+siteslist$Norm.D2Entropy)/2
siteslist$EntDiff<-siteslist$E_AvgEnt-siteslist$D_AvgEnt

condensedEntropies<-siteslist[1:2]
condensedEntropies$Gene_Symbol<-siteslist$Gene_Symbol
condensedEntropies$EntDiff<-siteslist$EntDiff
condensedEntropies$Euclidean<-siteslist$Euclidean
#condensedEntropies$Harmonic<-siteslist$Harmonic
colnames(condensedEntropies)[2]<-"NumSites"
ordered<-condensedEntropies[order(condensedEntropies$EntDiff,decreasing=TRUE),]

write_delim(data.frame(condensedEntropies$TransID, condensedEntropies$EntDiff), "refGene_tx_names.txt", delim="\t")

numtop<-sum(condensedEntropies$EntDiff>0)
top<-condensedEntropies[order(condensedEntropies$EntDiff,decreasing=TRUE),][1:numtop,]
numbot<-sum(condensedEntropies$EntDiff<0)
bot<-condensedEntropies[order(condensedEntropies$EntDiff,decreasing=FALSE),][1:numbot,]

stemgenes<-"Abcb5 Abcg2 Alcam Aldh1a1 Atm Atxn1 Axl Bmi1 Bmp7 Cd24 Cd34 Cd38 Cd44 Chek1 Ddr1 Dkk1 Dll1 Dll4 Dnmt1 Egf Eng Epcam Erbb2 Etfa Fgfr2 Flot2 Foxa2 Foxp1 Fzd7 Gata3 Gata4 Gsk3b Hdac1 Id1 Ikbkb Itga2 Itga4 Itga6 Itgb1 Jag1 Jak2 Kit Kitlg Klf17 Klf4 Lats1 Lin28a Lin28b Maml1 Mertk Ms4a1 Muc1 Myc Mycn Nanog Nfkb1 Nos2 Notch1 Notch2 Pecam1 Plat Plaur Pou5f1 Prom1 Ptch1 Ptprc Sav1 Sirt1 Smo Snai1 Snai2 Sox2 Stat3 Taz Tgfbr1 Thy1 Twist1 Twist2 Wee1 Wnt1 Wwc1 Yap1 Zeb1 Zeb2 Actb B2m Hprt1 Ldha Rplp1"
stemgenes<-data.frame(str_split(toupper(stemgenes), " "))
colnames(stemgenes)<-"Gene_Symbol"
stemgenesincommonIDs<-condensedEntropies$TransID[which(condensedEntropies$Gene_Symbol %in% stemgenes$Gene_Symbol)]
#stemgenesincommon<-intersect(stemgenes$Gene_Symbol, siteslist$Gene_Symbol)

#twist snail slug notch1 zeb2 ptch1 key drivers of emt-like transition
stemGeneRanks<-sort(which(ordered$TransID %in% stemgenesincommonIDs))
stemgenesEnt<-ordered[stemGeneRanks,]
stemgenesEnt$rank<-stemGeneRanks
#stemgenesEnt$Location<-allcompgb2$m6A_location[match(allcompgb2$Gene_Symbol, stemgenesincommon)]
