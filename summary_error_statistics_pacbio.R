\documentclass{article}

\begin{document}
\SweaveOpts{concordance=TRUE}

<<summary error statistics pacbio>>=
directory <-"/Users/alice/Desktop/projects/kinetics/errors_will/plot_summary_statistics"
setwd(directory)
par(mfrow=c(3,3))
filenames <- list.files(pattern = "*.mf", full.names = FALSE)
filenames <- filenames[grep("EmptyTmp", filenames, invert = TRUE)]
  
for (file in filenames) {
  print(file)
  emptyFile <- gsub(".mf", ".mfEmptyTmp", file)
  
  data <- read.table(
    file,
    header = F,
    col.names = c(
      "chr",
      "start",
      "end",
      "percTOTAL",
      "percMISMATCHES",
      "percINSERTION",
      "percDELETION",
      "TOTAL",
      "MISMATCHES",
      "INSERTION",
      "DELETION"
    ),
    sep = " ",
    fill = TRUE
  )
  
  dataEmpty <- read.table(
    emptyFile,
    header = F,
    col.names = c(
      "chr",
      "start",
      "end",
      "percTOTAL",
      "percMISMATCHES",
      "percINSERTION",
      "percDELETION",
      "TOTAL",
      "MISMATCHES",
      "INSERTION",
      "DELETION"
    ),
    sep = " ",
    fill = TRUE
  )
  
  #PRINT MOTIF FILE
  dfMotif <-
    cbind(
      sum(data$TOTAL, na.rm = TRUE),
      sum(data$MISMATCHES, na.rm = TRUE),
      sum(data$INSERTION, na.rm = TRUE),
      sum(data$DELETION, na.rm = TRUE),
      sum(data$MISMATCHES, na.rm = TRUE)/sum(data$TOTAL, na.rm = TRUE),
      sum(data$INSERTION, na.rm = TRUE)/sum(data$TOTAL, na.rm = TRUE),
      sum(data$DELETION, na.rm = TRUE)/sum(data$TOTAL, na.rm = TRUE)
    )
  
  #PRINT CONTROL FILE
  dfControl <-
    cbind(
      sum(dataEmpty$TOTAL, na.rm = TRUE),
      sum(dataEmpty$MISMATCHES, na.rm = TRUE),
      sum(dataEmpty$INSERTION, na.rm = TRUE),
      sum(dataEmpty$DELETION, na.rm = TRUE),
      sum(dataEmpty$MISMATCHES, na.rm = TRUE)/sum(dataEmpty$TOTAL, na.rm = TRUE),
      sum(dataEmpty$INSERTION, na.rm = TRUE)/sum(dataEmpty$TOTAL, na.rm = TRUE),
      sum(dataEmpty$DELETION, na.rm = TRUE)/sum(dataEmpty$TOTAL, na.rm = TRUE)
    )
  
  df <- data.matrix(rbind(dfMotif, dfControl))
  colnames(df) <- c("TOTAL",
                    "MISMATCHES",
                    "INSERTION",
                    "DELETION",
                    "fracMISMATCHES",
                    "fracINSERTION",
                    "fracDELETION")
  rownames(df) <- c("motifs","controls")
  print(df)
  barplot(df[,5:7],beside=TRUE,col=c("firebrick1","deepskyblue"),main=file)

}
@



\end{document}