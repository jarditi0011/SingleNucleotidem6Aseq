##completing annotations
## should probably be combined with the multicore script

features<-read.csv("forFeatureCalc.arraystar.csv")
symbolIDs<-data.frame(symbol=features$gene_id, id=features$transcript_id)
symbolIDs<-unique(symbolIDs)
#View(features)

siteslist_coding<-siteslist[str_detect(siteslist$TransID, "NM\\_"),]
allcompgb2_coding<-allcompgb2[str_detect(allcompgb2$TransID, "NM\\_"),]

##script to determine gene lengths for each tx in features

get_length<- function(txid) {
  mytx<-features[which(features$transcript_id==txid & features$ElemType=="exon"),]
  widths<-mytx$End-mytx$Start+1
  sum(widths)
}

#lengths<-allcompgb2[c(1,3,25)]
#lengths<-unique(lengths)
#lengths$my_length<-sapply(lengths$TransID, get_length)
#lengths$error<-lengths$my_length-lengths$transcript_Length

get_my_data<- function(txid) {
  my.tx.features<-features[which(features$transcript_id==txid),]
  my.tx.features$widths<-my.tx.features$End-my.tx.features$Start+1
  
  my.tx.exons<-my.tx.features[which(my.tx.features$ElemType=="exon"),]
  len<-sum(my.tx.exons$widths)
  
  CDSstart<-0
  if (my.tx.features$Strand[1]=="+") {
    start_codon_start<-min(my.tx.features$Start[which(my.tx.features$ElemType=="start_codon")])
    exonsIn5utr<-my.tx.exons[which(my.tx.exons$Start<start_codon_start),]
    if (nrow(exonsIn5utr)!=0) {
      CDSstart<-sum(exonsIn5utr$widths)
      if (max(exonsIn5utr$End)>start_codon_start) {
        CDSstart<-CDSstart-(max(exonsIn5utr$End)-start_codon_start+1)
      }
    }
  } else {
    start_codon_start<-max(my.tx.features$End[which(my.tx.features$ElemType=="start_codon")])
    exonsIn5utr<-my.tx.exons[which(my.tx.exons$End>start_codon_start),]
    if (nrow(exonsIn5utr)!=0) {
      CDSstart<-sum(exonsIn5utr$widths)
      if (min(exonsIn5utr$Start)<start_codon_start) {
        CDSstart<-CDSstart-(start_codon_start-min(exonsIn5utr$Start)+1)
      }
    }
  }
  
  my.tx.CDS<-my.tx.features[which(my.tx.features$ElemType=="CDS"),]
  CDSlen<-sum(my.tx.CDS$widths)
  
  return(c(CDSstart, CDSstart+CDSlen, len))
}

get_indices<-function(tx) {
  my.tx.indices<-strtoi(unlist(str_split(siteslist_coding$SiteLoci[which(siteslist_coding$TransID==tx)], pattern=", ")))
  return(my.tx.indices)
}

get_adjusted_indices<-function(tx) {
  
  data<-get_my_data(tx)
  cdsstart<-data[1]
  cdsend<-data[2]
  txend<-data[3]
  
  my.tx.indices<-get_indices(tx)
  numi<-length(my.tx.indices)
  for (num in 1:numi) {
    i<-my.tx.indices[num]
    if (i<cdsend){
      if (i<cdsstart) {
        i<-i/cdsstart
      } else {
        i<-(i-cdsstart)/(cdsend-cdsstart)+1
      }
    } else {
      i<-(i-cdsend)/(txend-cdsend)+2
    }
    my.tx.indices[num]<-i
  }
  return(my.tx.indices)
}

get_FCs<-function(tx) {
  my.tx.probs<-as.numeric(unlist(str_split(siteslist_coding$SiteProbs[which(siteslist_coding$TransID==tx)], pattern=", ")))
  return(my.tx.probs)
  }

plot_tx<-function(tx) {
  data<-get_my_data(tx)
  cdsstart<-data[1]
  cdsend<-data[2]
  txend<-data[3]
  ticks<-c(0, cdsstart, cdsend, txend)

  indices<-get_indices(tx)
  FCs<-get_FCs(tx)
  abundances<-data.frame(indices, FCs)
  
  gene<-siteslist_coding$Gene_Symbol[which(siteslist_coding$TransID==tx)]
  my_title<-paste(gene, ", ", as.character(txend), " bp", sep="")
  
    p<-ggplot(abundances, aes(indices, FCs, ymax = FCs, ymin = 0)) +
    geom_point(size=0.1) +
    geom_pointrange() +
    ylab("Log2 Fold Change") +
    theme(panel.background = element_blank(), 
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y=element_text(size=8),
          axis.ticks = element_blank()) +
    xlim(0, txend) +
    ylim(min(allcompgb2$Foldchange), max(allcompgb2$Foldchange)) +
    geom_hline(yintercept = 0) +
    geom_rect(mapping=aes(xmin=cdsstart, xmax=cdsend, ymin=max(allcompgb2$Foldchange)*-0.05, ymax=max(allcompgb2$Foldchange)*0.05)) +
    geom_text(aes(x=cdsstart+100, y=max(allcompgb2$Foldchange)*0.15, label="CDS")) +
    geom_text(aes(x=txend, y=max(allcompgb2$Foldchange)*0.1, label=gene, fontface="italic")) +
    geom_text(aes(x=txend, y=max(allcompgb2$Foldchange)*-0.1, label=paste(txend, "bp"))) 

  #ggsave("test.jpeg", 
  #       plot=p,
  #       device ="jpeg", 
  #       width=8, 
  #       height=2, 
  #       units = "in")
  
  ggsave(paste(gene, ".jpeg", sep=""), 
         plot=p,
         device ="jpeg", 
         width=8, 
         height=2, 
         units = "in")
}



plot_txs<-function(txs, title) {
  prevL<-length(txs)
  txs<-txs[str_detect(txs, "NM\\_")]
  newL<-length(txs)
  print(paste("Reduced from", prevL, "transcripts to", newL, "coding transcripts."))
  indices<-c()
  FCs<-c()
  ids<-c()
  
  for (num in 1:length(txs)){
    tx<-txs[num]
    my.indices<-get_adjusted_indices(tx)
    my.FCs<-get_FCs(tx)
    my.ids<-rep(tx, times=length(my.indices))
    ids<-c(ids, my.ids)
    indices<-c(indices, my.indices)
    FCs<-c(FCs, my.FCs)
  }
  
  ticks<-c(0, 1, 2, 3)
  abundances<-data.frame(ids, indices, FCs)
  
  my_title<-if (nchar(title)>10) {
    ""
  } else {
    title
  }
  
  p<-ggplot(abundances, aes(indices, FCs, ymax = FCs, ymin = 0)) +
    geom_point(size=0.1) +
    geom_pointrange() +
    ylab("Log2 Fold Change") +
    theme(panel.background = element_blank(), 
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y=element_text(size=8),
          axis.ticks = element_blank()) +
    xlim(0, 3) +
    ylim(min(allcompgb2$Foldchange), max(allcompgb2$Foldchange)) +
    geom_hline(yintercept = 0) +
    geom_rect(mapping=aes(xmin=1, xmax=2, ymin=max(allcompgb2$Foldchange)*-0.05, ymax=max(allcompgb2$Foldchange)*0.05)) +
    geom_text(aes(x=1.1, y=max(allcompgb2$Foldchange)*0.15, label="CDS")) +
    geom_text(aes(x=3, y=max(allcompgb2$Foldchange)*0.1, label=my_title, fontface="italic")) 
  
  ggsave(paste(title, "txs.jpeg", sep=""), 
         plot=p,
         device ="jpeg", 
         width=8, 
         height=2, 
         units = "in")
}

plot_pie<-function(txs, title) {
  
  txs<-remove_nc(txs) 
  
  locations<-allcompgb2_coding$m6A_location[which(allcompgb2_coding$TransID %in% txs)]
  loc<-data.frame(table(locations))
  
  # Compute the position of labels
  loc <- loc %>% 
    arrange(desc(locations)) %>%
    mutate(prop = Freq / sum(loc$Freq) *100) %>%
    mutate(ypos = cumsum(prop)- 0.5*prop )
  
  # Basic piechart
  pie<- ggplot(loc, aes(x="", y=prop, fill=locations)) +
    geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0) +
    theme_void() + 
    theme(legend.position="none") +
    geom_text(aes(y = ypos, label = paste(locations, "\n", Freq)), color = "white", size=3) +
    scale_fill_brewer(palette="Set1") +
    ggtitle(title)
  
  ggsave(paste(title, "pie.jpeg"), 
         plot=pie,
         device ="jpeg", 
         width=4, 
         height=4, 
         units = "in")
}

remove_nc<-function(txs) {
  prevL<-length(txs)
  txs<-txs[str_detect(txs, "NM\\_")]
  newL<-length(txs)
  print(paste("Reduced from", prevL, "transcripts to", newL, "coding transcripts."))
  return(txs)
}

plot_tx_dist<- function(txs, title, boxscale, cdsscale) {
  
  txs<-remove_nc(txs)
  newL<-length(txs)
  
  plot_pie(txs, title)
  
  indices<-c()
  FCs<-c()
  ids<-c()
  
  for (num in 1:length(txs)){
    tx<-txs[num]
    my.indices<-get_adjusted_indices(tx)
    my.FCs<-get_FCs(tx)
    my.ids<-rep(tx, times=length(my.indices))
    ids<-c(ids, my.ids)
    indices<-c(indices, my.indices)
    FCs<-c(FCs, my.FCs)
  }
  
  ticks<-c(0, 1, 2, 3)
  abundances<-data.frame(ids, indices, FCs)
  
  my_title<-title
  
  p2<-ggplot(abundances, aes(indices)) +
    geom_freqpoly(mapping=aes(x=indices), 
                   inherit.aes = FALSE, 
                   color="black",
                   bins=floor(newL*(300/4346))) +
    ylab("Number of Sites") +
    xlim(0,3) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y=element_text(size=8),
          axis.ticks = element_blank()) +
    geom_segment(aes(x=0,xend=3,y=0,yend=0)) +
    geom_text(aes(x=1.075, y=cdsscale, label="CDS")) +
    ggtitle(title) +
    geom_rect(mapping=aes(xmin=1, xmax=2, ymin=-boxscale, ymax=boxscale))
  
  ggsave(paste(title, "dist.jpeg", sep=""), 
         plot=p2,
         device ="jpeg", 
         width=8, 
         height=4, 
         units = "in")
  
  ggsave("test.jpeg", 
         plot=p2,
         device ="jpeg", 
         width=8, 
         height=4, 
         units = "in")
}

get_adjusted_index<-function(tx, i) {
  
  data<-get_my_data(tx)
  cdsstart<-data[1]
  cdsend<-data[2]
  txend<-data[3]
  
  if (i<cdsend){
    if (i<cdsstart) {
      newi<-i/cdsstart
    } else {
      newi<-(i-cdsstart)/(cdsend-cdsstart)+1
    }
  } else {
    newi<-(i-cdsend)/(txend-cdsend)+2
  }
return(newi)
}

#mysites is a dataframe with the same columns as allcompgb2
plot_txs_bysite<- function(mysites, title, boxscale, cdsscale) {
  #mysites<-allcompgb2_coding[which(allcompgb2_coding$TransID %in% txs),]
  
  abundances<-data.frame(ids=mysites$TransID, 
                         prep_indices=mysites$m6A_transcript_location, 
                         FCs=mysites$Foldchange)
  abundances$indices<-mapply(get_adjusted_index, 
                             abundances$ids, 
                             abundances$prep_indices)
  newL<-nrow(abundances)

  ticks<-c(0, 1, 2, 3)
  my_title<-title
  
  p2<-ggplot(abundances, aes(indices)) +
    geom_freqpoly(mapping=aes(x=indices), 
                  inherit.aes = FALSE, 
                  color="black",
                  bins=floor(newL*(300/4346))) +
    ylab("Number of Sites") +
    xlim(0,3) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y=element_text(size=8),
          axis.ticks = element_blank()) +
    geom_segment(aes(x=0,xend=3,y=0,yend=0)) +
    geom_text(aes(x=1.075, y=cdsscale, label="CDS")) +
    ggtitle(title) +
    geom_rect(mapping=aes(xmin=1, xmax=2, ymin=-boxscale, ymax=boxscale))
  
  ggsave(paste(title, "dist.jpeg", sep=""), 
         plot=p2,
         device ="jpeg", 
         width=8, 
         height=4, 
         units = "in")
  
  locations<-mysites$m6A_location
  loc<-data.frame(table(locations))
  
  # Compute the position of labels
  loc <- loc %>% 
    arrange(desc(locations)) %>%
    mutate(prop = Freq / sum(loc$Freq) *100) %>%
    mutate(ypos = cumsum(prop)- 0.5*prop )
  
  # Basic piechart
  pie<- ggplot(loc, aes(x="", y=prop, fill=locations)) +
    geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0) +
    theme_void() + 
    theme(legend.position="none") +
    geom_text(aes(y = ypos, label = paste(locations, "\n", Freq)), color = "white", size=3) +
    scale_fill_brewer(palette="Set1") +
    ggtitle(title)
  
  ggsave(paste(title, "pie.jpeg"), 
         plot=pie,
         device ="jpeg", 
         width=4, 
         height=4, 
         units = "in")
  
  my_title<-if (nchar(title)>10) {
    ""
  } else {
    title
  }
  
  p<-ggplot(abundances, aes(indices, FCs, ymax = FCs, ymin = 0)) +
    geom_point(size=0.1) +
    geom_pointrange() +
    ylab("Log2 Fold Change") +
    theme(panel.background = element_blank(), 
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y=element_text(size=8),
          axis.ticks = element_blank()) +
    xlim(0, 3) +
    ylim(min(allcompgb2$Foldchange), max(allcompgb2$Foldchange)) +
    geom_hline(yintercept = 0) +
    geom_rect(mapping=aes(xmin=1, xmax=2, ymin=max(allcompgb2$Foldchange)*-0.05, ymax=max(allcompgb2$Foldchange)*0.05)) +
    geom_text(aes(x=1.1, y=max(allcompgb2$Foldchange)*0.15, label="CDS")) +
    geom_text(aes(x=3, y=max(allcompgb2$Foldchange)*0.1, label=my_title, fontface="italic")) 
  
  ggsave(paste(title, "txs.jpeg", sep=""), 
         plot=p,
         device ="jpeg", 
         width=8, 
         height=2, 
         units = "in")
  
}

