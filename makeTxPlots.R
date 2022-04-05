##intense workflow
sitesordered<-allcompgb2[order(allcompgb2$Foldchange),]
sitesordered<-sitesordered[str_detect(sitesordered$TransID, "NM\\_"),]
top800bysite<-sitesordered[1:800,]
bot800bysite<-sitesordered[7860:8659,]
mid1600bysite<-sitesordered[3530:5129,]
ordered<-ordered[str_detect(ordered$TransID, "NM\\_"),]
top300<-ordered$TransID[1:300]
bot300<-ordered$TransID[(nrow(ordered)-299):nrow(ordered)]
mid600<-ordered$TransID[((nrow(ordered)/2)-299):((nrow(ordered)/2)+300)]

test<-read.csv("kinase binding  .csv", skip = 1, header=F)$V1

useme<-top800bysite
title<-"Top 800 sites by FC"

#plot_pie(useme, title)
plot_txs_bysite(useme, title, 1, 2.5)

plot_tx_dist(useme, title, 2,5)
plot_txs(useme, title)


stemtxs<-stemgenesEnt$TransID
numtxs<-length(stemtxs)
for (i in 1:numtxs) {
  plot_tx(stemtxs[i])
}
