##file containing plots of data

h<-hist(siteslist$EntDiff,
        main = "Distribution of Entropy Changes", 
        xlab = "Entropy Difference (Entacapone - DMSO)",
        #xaxp = c(-0.04, 0.3, 17),
        ylab = "Number of Transcripts")
abline(v=stemgenesEnt$EntDiff, col="red")
text(stemgenesEnt$EntDiff[1:5]+0.0075, c(2000, 2000, 2000, 1500, 1000), stemgenesEnt$Gene_Symbol[1:5], cex=0.6, srt=90)
text(0.2, 1250, "Stemness genes \n all increased \n entropy after \n treatment", col="red", cex=0.6)
#text(-0.2, 2000, paste("Other stemness genes that increased entropy: \n", str(stemgenesEnt$Gene_Symbol[6:14])), cex=0.5)
d<-density(siteslist$EntDiff)
lines(d,
      main = "Density of Entropy Changes", 
      xlab = "Entropy Difference (Entacapone - DMSO)", 
      ylab = "Density (~Number of Transcripts)")

plot(condensedEntropies$EntDiff, condensedEntropies$Euclidean,
     main="Entropy Difference vs L1 Norm", 
     xlab="Entropy Difference", 
     ylab="L1 Norm")
abline(lm(condensedEntropies$Mean ~ condensedEntropies$EntDiff))


topEuc<-condensedEntropies$Euclidean[which(condensedEntropies$EntDiff>0)]
topEnt<-condensedEntropies$EntDiff[which(condensedEntropies$EntDiff>0)]
title<-"Signed Euclidean Distance and Entropy Difference"

p<-ggplot(condensedEntropies, 
       aes(x=EntDiff, y=Euclidean)) + 
  theme(aspect.ratio=1, 
        plot.title = element_text(size = 11, 
                                  hjust = 0.5),
        axis.title.x=element_text(size=8),
        axis.title.y=element_text(size=8),) +
  geom_point() +
  #ggtitle("\n Signed Euclidean Distance \n and Entropy Difference") + 
  xlab("Entropy Difference") +
  ylab("Signed Euclidean Distance") +
  geom_smooth(method=lm , color="gray22", se=TRUE) +
  stat_poly_eq(aes(x=EntDiff, y=Euclidean),#aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE, 
               na.rm=TRUE, 
               coef.keep.zeros = TRUE) +
  geom_point(data=ordered[1:300,],
             aes(x=EntDiff,y=Euclidean),
             color='blue', 
             size=0.5) +
  geom_point(data=ordered[(nrow(ordered)-299):nrow(ordered),],
             aes(x=EntDiff,y=Euclidean),
             color='red', 
             size=0.5) +
  coord_fixed()
p<-ggExtra::ggMarginal(p, size=10)
plot(p)

ggsave(paste(title, "300.pdf"),
       plot=p,
       height = 5,
       width=5,
       units="in")

##########STEM

coolstemGenes<-c("NOTCH1", "GSK3B", "DLL1", "ALCAM", "PTCH1", "ZEB2", "MAML1", "ATXN1", "HDAC1")

p<-ggplot(condensedEntropies, 
          aes(x=EntDiff, y=Euclidean)) + 
  theme(aspect.ratio=1, 
        plot.title = element_text(size = 11, 
                                  hjust = 0.5),
        axis.title.x=element_text(size=8),
        axis.title.y=element_text(size=8),) +
  geom_point() +
  #ggtitle("\n Signed Euclidean Distance \n and Entropy Difference") + 
  xlab("Entropy Difference") +
  ylab("Signed Euclidean Distance") +
  geom_smooth(method=lm , color="gray22", se=TRUE) +
  stat_poly_eq(aes(x=EntDiff, y=Euclidean),#aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE, 
               na.rm=TRUE, 
               coef.keep.zeros = TRUE) +
  geom_point(data=stemgenesEnt, 
             aes(x=EntDiff,y=Euclidean), 
             color='tomato4', 
             size=0.5) +
  geom_label_repel(aes(label=ifelse(
    (Gene_Symbol %in% coolstemGenes),
    as.character(Gene_Symbol),'')),
    size=2.5,
    #label.size = NA,
    max.overlaps = Inf,
    color="tomato4") +
  coord_fixed()
p<-ggExtra::ggMarginal(p, size=10)
#plot(p)

ggsave(paste(title, "stem.pdf"),
       plot=p,
       height = 5,
       width=5,
       units="in")

########### complexes
repc<-c("EP300", "MAML1")
actc<-c("HDAC1", "NCOR1")
repdf<-stemgenesEnt[which(stemgenesEnt$Gene_Symbol %in% repc),]
actdf<-stemgenesEnt[which(stemgenesEnt$Gene_Symbol %in% actc),]


p<-ggplot(condensedEntropies, 
          aes(x=EntDiff, y=Euclidean)) + 
  theme(aspect.ratio=1, 
        plot.title = element_text(size = 11, 
                                  hjust = 0.5),
        axis.title.x=element_text(size=8),
        axis.title.y=element_text(size=8),) +
  geom_point() +
  #ggtitle("\n Signed Euclidean Distance \n and Entropy Difference") + 
  xlab("Entropy Difference") +
  ylab("Signed Euclidean Distance") +
  geom_smooth(method=lm , color="gray22", se=TRUE) +
  stat_poly_eq(aes(x=EntDiff, y=Euclidean),#aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE, 
               na.rm=TRUE, 
               coef.keep.zeros = TRUE) +
  geom_point(data=repdf, 
             aes(x=EntDiff,y=Euclidean), 
             color='red4', 
             size=0.5) +
  geom_point(data=actdf, 
             aes(x=EntDiff,y=Euclidean), 
             color='blue4', 
             size=0.5) +
  geom_label_repel(aes(label=ifelse(
    (Gene_Symbol %in% repc),
    as.character(Gene_Symbol),'')),
    size=2.5,
    #label.size = NA,
    max.overlaps = Inf,
    color="red4") +
  geom_label_repel(aes(label=ifelse(
    (Gene_Symbol %in% actc),
    as.character(Gene_Symbol),'')),
    size=2.5,
    #label.size = NA,
    max.overlaps = Inf,
    color="blue4") +
  coord_fixed()
p<-ggExtra::ggMarginal(p, size=10)
#plot(p)

ggsave(paste(title, "complex.pdf"),
       plot=p,
       height = 5,
       width=5,
       units="in")

###"Entropy as a function of probability"
title<-"Entropy as a function of probability"

entropy<-function(x) {
  if (x==0 | x==1) {
    0
  } else {
  -1*(x*log(x)+(1-x)*log(1-x))
  }
}
vals<-seq(0,1,0.1)

ents<-data.frame(site_probability=vals, site_entropy=sapply(vals, entropy))

p<-ggplot(ents, aes(x=site_probability, y=site_entropy)) +
  geom_point() +
  geom_line() +
  ggtitle(title) +
  xlab("Site Probability") +
  ylab("Site Entropy") +
  theme(plot.title = element_text(size = 11, 
                                  hjust = 0.5),
        axis.title.x=element_text(size=8),
        axis.title.y=element_text(size=8),) +
  coord_fixed()

ggsave(paste(title, ".pdf"),
       plot=p,
       height = 5,
       width=5,
       units="in", 
       device = "pdf")

#####

