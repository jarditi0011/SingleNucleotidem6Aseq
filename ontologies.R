

top300<-as.list(read_delim("C:/Users/jardi/OneDrive/Documents/TapinosLab/singlenucm6a/top300.txt", delim="\n"))$top300.Gene_Symbol
#bot300<-as.list(read_delim("bot300.txt", delim="\n"))$bot300.Gene_Symbol

#setEnrichrSite("Enrichr")
#websiteLive <- TRUE
#dbs <- listEnrichrDbs()
#if (is.null(dbs)) websiteLive <- FALSE
#if (websiteLive) head(dbs)
####COMMENT ALL ABOVE OUT AFTER 1ST RUN

write_ontology<-function(rank) {
  write_csv(data.frame(str_split(df$Genes[rank], pattern=";")), paste(df$shortName[rank],toporbot,".csv"))
}

useTop<-FALSE
databasenum=4
showTerms = 20 
numChar = 40 
y = "Count" 
orderBy = "P.value"

mydbs <- c("GO_Molecular_Function_2021", 
           "GO_Cellular_Component_2021", 
           "GO_Biological_Process_2021", 
           "KEGG_2021_Human")
dbsnames<-c("Molecular function of", 
            "Cellular components enriched by", 
            "Biological processes involving the", 
            "KEGG pathways involving the")

if (useTop) {
  dataToUse<-top300
  toporbot<-"top"
} else {
  dataToUse<-bot300
  toporbot<-"bottom"
}


## Given a Enrichr output, order and subset criteria, returns a data frame accordingly
.enrichment_prep_df <- function(df, showTerms, orderBy) {
  
  if(is.null(showTerms)) {
    showTerms = nrow(df)
  } else if(!is.numeric(showTerms)) {
    stop(paste0("showTerms '", showTerms, "' is invalid."))
  }
  
  Annotated <- as.numeric(sub("^\\d+/", "", as.character(df$Overlap)))
  Significant <- as.numeric(sub("/\\d+$", "", as.character(df$Overlap)))
  
  # Build data frame
  df <- cbind(df, data.frame(Annotated = Annotated, Significant = Significant,
                             stringsAsFactors = FALSE))
  
  # Order data frame (P.value or Combined.Score)
  if(orderBy == "Combined.Score") {
    idx <- order(df$Combined.Score, decreasing = TRUE)
  } else {
    idx <- order(df$P.value, decreasing = FALSE)
  }
  df <- df[idx,]
  
  # Subset to selected number of terms
  if(showTerms <= nrow(df)) {
    df <- df[1:showTerms,]
  }
  
  return(df)
}

if (websiteLive) {
  enriched <- enrichr(dataToUse, mydbs)
}
#if (websiteLive) enriched[["GO_Biological_Process_2021"]]
#if (websiteLive) plotEnrich(enriched[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")

df=enriched[[databasenum]] 
df <- .enrichment_prep_df(df, showTerms, orderBy)

if (databasenum!=4){
  shortName<-strsplit(df$Term, split="\\(GO")
  shortName<-unlist(shortName)[ c(TRUE,FALSE) ]
} else {
  shortName<-df$Term
}
#shortName <- paste(substr(df$Term, 1, numChar), ifelse(nchar(df$Term) > 
#numChar, "...", ""), sep = "")

df$shortName = shortName#df$Term #use df$Term if trying to access KEGG pathways
df$shortName <- factor(df$shortName, levels = rev(unique(df$shortName)))
df$Ratio <- df$Significant/df$Annotated

if (orderBy == "Combined.Score") {
  fill <- "Combined.Score"
} else {
  fill <- "P.value"
}
if (y != "Ratio") {
  y <- "Significant"
}

map <- aes_string(x = "shortName", y = y, fill = fill)
xlab <- "Enriched terms"
if (y == "Ratio") {
  ylab <- "Gene ratio"
} else {
  ylab <- "Gene count"
}
title <- paste(dbsnames[databasenum], toporbot, "300 transcripts")

p <- ggplot(df, map) + geom_bar(stat = "identity") + 
  coord_flip() + theme_bw()
if (orderBy == "Combined.Score") {
  p <- p + scale_fill_continuous(low = "blue", high = "red") + 
    guides(fill = guide_colorbar(title = "Combined Score", 
                                 reverse = FALSE))
} else {
  p <- p + scale_fill_continuous(low = "red", high = "blue") + 
    guides(fill = guide_colorbar(title = "P value", 
                                 reverse = TRUE))
}
p <- p + theme(axis.text.x = element_text(colour = "black", 
                                          vjust = 1), axis.text.y = element_text(colour = "black", 
                                                                                 hjust = 1), axis.title = element_text(color = "black", 
                                                                                                                       margin = margin(10, 5, 0, 0)), axis.title.y = element_text(angle = 90))
p <- p + xlab(xlab) + ylab(ylab) + ggtitle(title)
plot(p)
View(df)



#ggsave("test.jpeg", 
#       plot=p, 
#       device ="jpeg", 
#       width=8, 
#       height=4, 
#       units = "in")