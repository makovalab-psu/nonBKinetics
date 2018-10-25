#/galaxy/home/biomonika/R-3.2.4revised/bin/Rscript generateErrorStatistics.R motifFile pathTocmpFile

library("h5r", lib.loc = "/galaxy/home/biomonika/R-3.2.4revised_nn/library/")
library("pbh5", lib.loc = "/galaxy/home/biomonika/R-3.2.4revised_nn/lib")
library("data.table", lib.loc = "/galaxy/home/biomonika/R-3.2.4revised_nn/lib")
options(show.error.locations = TRUE)

args <- commandArgs(trailingOnly = TRUE)
pathTocmpFile = args[1]
reference = args[2]
motif_folder = args[3]
out_location = args[4]

setwd(out_location)

cmpH5file = paste(pathTocmpFile, "chr", reference, "_P6.cmp.h5", sep = "")
cmp = PacBioCmpH5(cmpH5file)

getWindow <- function(arguments) {
  #print(arguments)
  chr <- arguments[1]
  start <- as.numeric(arguments[2])
  end <- as.numeric(arguments[3])
  
  rates <- (getErrorRate(start, end))
  percErrorTotal <- (rates$mismatchRows + rates$insertionRows + rates$deletionRows) / rates$totalRows
  percErrorIns <- rates$insertionRows / rates$totalRows
  percErrorDel <- rates$deletionRows / rates$totalRows
  percErrorMism <- rates$mismatchRows / rates$totalRows
  
  return(
    paste(
      reference,
      start,
      end,
      percErrorTotal,
      percErrorMism,
      percErrorIns,
      percErrorDel,
      rates$totalRows,
      rates$mismatchRows,
      rates$insertionRows,
      rates$deletionRows,
      sep = " "
    )
  )
  #rows<-rbind(rows,t)
  #write.table(t, file=filename, col.names = FALSE,row.names = FALSE,quote=FALSE,sep="\t",append=TRUE)
  
}

getErrorRate <- function(start, end) {
  # start
  # end
  subreads <- getReadsInRange(cmp, paste0("chr",reference), start, end)
  
  if (length(subreads)>0) {
    w <- getByTemplatePosition(cmp, idx = subreads)
    w <- as.data.table(w)
    w <- w[position >= start & position <= end]
    
    if (nrow(w) > 0) {
      totalRows <- nrow(w)
      
      
    insertion_lengths <-
      unlist(lapply(unique(w$idx), function(x)
        (rle(w[idx == x]$ref)$lengths[rle(w[idx == x]$ref)$values == "-"])))
    if (!length(insertion_lengths) > 0) {
      insertion_lengths <- NA
    }
    
    deletion_lengths <-
      unlist(lapply(unique(w$idx), function(x)
        (rle(w[idx == x]$read)$lengths[rle(w[idx == x]$read)$values == "-"])))
    if (!length(deletion_lengths) > 0) {
      deletion_lengths <- NA
    }
      
      insertionRows <- length(insertion_lengths)
      deletionRows <- length(deletion_lengths)
      mismatchRows <-
        nrow(w[as.character(w$ref) != as.character(w$read) &
                 (as.character(w$ref) != "-") & (as.character(w$read) != "-"), ])
      
      return(list(
        totalRows=totalRows,
        mismatchRows=mismatchRows,
        insertionRows=insertionRows,
        deletionRows=deletionRows
      ))
    }
  }
}

processMotif <- function(motif) {
  print("ITERATION")
  print(reference)
  print(motif)
  
  filename <-
    paste(reference,
          ".ERRORS",
          basename(motif_folder),
          ".",
          basename(motif),
          ".txt",
          sep = "")
  motifFile <- paste(motif_folder, "/", motif, sep = "")
  
  if (file.exists(motifFile)) {
    #.mf file exists
    coordinates <-
      read.table(motifFile,col.names = paste0("V",seq_len(max(count.fields(motifFile)))), fill = TRUE)[, 1:3] #read only first three columns
    coordinates <- as.data.table(coordinates)
    setkey(coordinates, "V1", "V2", "V3")
    print(dim(coordinates))
    coordinates <-
      subset(coordinates, V1 == reference) #subset only to specific chromosome #paste("chr",reference,sep="")
    print(dim(coordinates))
    
    #rows<-NULL
    rows <- apply(coordinates, 1, function(x)
      getWindow(x))

    if (file.exists(filename)) {
      file.remove(filename) #file with results already exists, remove before writing
    }

    system.time(
      write.table(
        as.matrix(rows),
        file = filename,
        col.names = FALSE,
        row.names = FALSE,
        quote = FALSE,
        sep = "\t",
        append = TRUE
      )
    )
  } else {
    print("File does not exist, skipping.")
  }
}

list<-list.files(path=motif_folder,pattern = "*.mf", full.names = FALSE)

print("NON-PARALLEL")
ptm <- proc.time()
for (motif in list) { 
  processMotif(motif) #PROCESS EACH MOTIF INDIVIDUALLY
}
proc.time() - ptm #stop timer

print("Done.")