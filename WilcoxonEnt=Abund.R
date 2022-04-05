#this confirms that difference in entropy indeed distinguishes site abundances

topForWilcox<-allcompgb2[which(allcompgb2$Gene_Symbol %in% top$Gene_Symbol[1:500]),]
TopEForWilcox<-sapply(sapply(topForWilcox[11], as.character), as.numeric)
TopDForWilcox<-sapply(sapply(topForWilcox[12], as.character), as.numeric)
wilcox.test(TopEForWilcox, TopDForWilcox)
median(topForWilcox$E.Mean.pm6A)
median(topForWilcox$D.Mean.pm6A)
#class(colnames(topForWilcox[11:12]))

BotForWilcox<-allcompgb2[which(allcompgb2$Gene_Symbol %in% bot$Gene_Symbol[1:500]),]
BotEForWilcox<-sapply(sapply(BotForWilcox[11], as.character), as.numeric)
BotDForWilcox<-sapply(sapply(BotForWilcox[12], as.character), as.numeric)
wilcox.test(BotEForWilcox, BotDForWilcox)
median(BotForWilcox$E.Mean.pm6A)
median(BotForWilcox$D.Mean.pm6A) 
#TODO plot in ggplot, boxplot, violin plot

