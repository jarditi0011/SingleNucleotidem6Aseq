##annotating multicore
library(foreach)
library(doParallel)
#C:/Users/jardi/OneDrive/Documents/TapinosLab/Annotataed HG/hg19.refGene.gtf"
genome<-read.table("hg19.refGene2014.gtf", sep="\t")
#realgenome<-genome
#genome<-realgenome[1:1000,]
allcompgb2<-read.csv("GB2ED.csv")
colnames(allcompgb2)[which(names(allcompgb2)=="ï..TransID")]<-"TransID"
keys<-unique(allcompgb2$TransID)
#genome<-realgenome
colnames(genome)<-c("Chr", "Src", "ElemType", "Start", "End", "Score", "Strand", "Phase", "Info")

#genome$gene_id<-""
#genome$transcript_id<-""
#genome$exon_number<-""
#genome$exon_id<-""
#genome$gene_name<-""

registerDoParallel(cores = detectCores()-1)
cat(paste("Hyper-thread registered:", getDoParRegistered(), 
          "\n"))
cat(paste("Using", getDoParWorkers(), 
          "thread(s) to count unpack annotations...\n"))
elems<-nrow(genome)
start_time<-Sys.time()
print(paste("Started at ", starttime))

info<-foreach(i = 1:elems) %dopar% {
  val<-genome[i,9]
  splitstr<-strsplit(val, split="; ")
  elemslist<-unlist(sapply(splitstr, function(x) strsplit(x, split=" ")))
  
  gene_id<-elemslist[match("gene_id", elemslist)+1]
  transcript_id<-elemslist[match("transcript_id", elemslist)+1]
  exon_number<-elemslist[match("exon_number", elemslist)+1]
  exon_id<-elemslist[match("exon_id", elemslist)+1]
  
  gene_name<-gsub(";", "", elemslist[match("gene_name", elemslist)+1])
  genome[i, 10:14]<-c(gene_id, transcript_id, exon_number, exon_id, gene_name)
}
end_time <- Sys.time()
print(paste("Time used to unpack:", 
            difftime(end_time, start_time, 
                     units = "mins"), "mins"))
info <- as.data.frame(do.call(rbind, info))
colnames(info)<-c("gene_id","transcript_id","exon_number","exon_id","gene_name")
genome<-cbind(genome, info)

cat("getting relevant features... \n")
start_time<-Sys.time()
mytxs<-genome[which(genome$transcript_id %in% keys),]
end_time<-Sys.time()
print(paste("Time used to filter:", 
            difftime(end_time, start_time, 
                     units = "mins"), "mins"))

mytxsrel<-mytxs[c(1,3,4,5,7,10,11,12,13,14)]

write.csv(genome, "hg19.refGene2014.unpacked.gtf")
write.csv(mytxsrel, "forFeatureCalc2014.csv")
