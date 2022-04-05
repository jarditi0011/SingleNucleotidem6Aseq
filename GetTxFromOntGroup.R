path<-"C:/Users/jardi/OneDrive/Documents/TapinosLab/singlenucm6a/OntologyGroups/"
files <- list.files(path=path, pattern="*.csv", full.names=TRUE, recursive=FALSE)
lapply(files, function(file) {
  
  genes<-read.csv(file, skip = 1, header=F)$V1
  dir<-str_remove(file, ".csv")
  title<-str_remove(dir, path)
  genes<-symbols2ids(genes)
  genes<-genes[str_detect(genes, "NM\\_")]
  numgenes<-length(genes)
  dir.create(dir)
  setwd(dir)
  print(paste("working on", title))
  for (num in 1:numgenes) {
    gene<-genes[num]
    #print(gene)
    plot_tx(gene)
  }
  print("complete!")
  setwd(path)
})
