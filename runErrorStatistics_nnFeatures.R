#/galaxy/home/biomonika/R-3.2.4revised/bin/Rscript generateErrorStatistics.R motifFile pathTocmpFile

library("h5r", lib.loc = "/galaxy/home/biomonika/R-3.2.4revised_nn/library/")
library("pbh5", lib.loc = "/galaxy/home/biomonika/R-3.2.4revised_nn/lib")
library("data.table", lib.loc = "/galaxy/home/biomonika/R-3.2.4revised_nn/lib")
options(show.error.locations = TRUE)

setwd("/nfs/brubeck.bx.psu.edu/scratch6/monika/pac_errors/latest")
args <- commandArgs(trailingOnly = TRUE)
pathTocmpFile = args[1]
reference = args[2]
folder = args[3]

cmpH5file = paste(pathTocmpFile, "chr", reference, "_P6.cmp.h5", sep = "")
cmp = PacBioCmpH5(cmpH5file)

##get IPDs by template position
res = getByTemplatePosition(cmp, f = getIPD)
RES <- as.data.table(res)
setkey(RES, "position", "read", "ref")


getWindow <- function(arguments) {
  #print(arguments)
  chr <- arguments[1]
  start <- as.numeric(arguments[2])
  end <- as.numeric(arguments[3])
  
  rates <- (getErrorRate(start, end))
  percErrorTotal <- (rates[[2]] + rates[[3]] + rates[[4]]) / rates[[1]] *
    100
  percErrorIns <- (rates[[2]]) / rates[[1]] * 100
  percErrorDel <- (rates[[3]]) / rates[[1]] * 100
  percErrorMism <- (rates[[4]]) / rates[[1]] * 100
  
  return(
    paste(
      reference,
      start,
      end,
      rates[[1]],
      rates[[2]],
      rates[[3]],
      rates[[4]],
      percErrorTotal,
      percErrorIns,
      percErrorDel,
      percErrorMism,
      paste0("$", rates[[5]], "$", rates[[6]]),
      sep = " "
    )
  )
  #rows<-rbind(rows,t)
  #write.table(t, file=filename, col.names = FALSE,row.names = FALSE,quote=FALSE,sep="\t",append=TRUE)
  
}

getErrorRate <- function(start, end) {
  # start
  # end
  w <-
    RES[position >= start &
          position < end] #cut out relevant portion from the data
  dim(w)
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
    
    insertionRows <- nrow(w[as.character(w$ref) == "-", ])
    deletionRows <- nrow(w[as.character(w$read) == "-", ])
    mismatchRows <-
      nrow(w[as.character(w$ref) != as.character(w$read) &
               (as.character(w$ref) != "-") & (as.character(w$read) != "-"), ])
    
    return(list(
      totalRows,
      insertionRows,
      deletionRows,
      mismatchRows,
      list(as.numeric(insertion_lengths)),
      list(as.numeric(deletion_lengths))
    ))
  }
}

processMotif <- function(motif) {
  print("ITERATION")
  print(reference)
  print(motif)
  
  filename <-
    paste(reference,
          ".ERRORS",
          folder,
          ".",
          basename(motif),
          ".txt",
          sep = "")
  motifFile <- paste(folder, "/", motif, sep = "")
  
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

listMotifs <-
  c(
    "AAACnFeatureOnly.mf",
    "AAAGnFeatureOnly.mf",
    "AAATnFeatureOnly.mf",
    "AACCnFeatureOnly.mf",
    "AACGnFeatureOnly.mf",
    "AACnFeatureOnly.mf",
    "AACTnFeatureOnly.mf",
    "AAGCnFeatureOnly.mf",
    "AAGGnFeatureOnly.mf",
    "AAGnFeatureOnly.mf",
    "AAGTnFeatureOnly.mf",
    "AATCnFeatureOnly.mf",
    "AATGnFeatureOnly.mf",
    "AATnFeatureOnly.mf",
    "AATTnFeatureOnly.mf",
    "ACAGnFeatureOnly.mf",
    "ACATnFeatureOnly.mf",
    "ACCCnFeatureOnly.mf",
    "ACCGnFeatureOnly.mf",
    "ACCnFeatureOnly.mf",
    "ACCTnFeatureOnly.mf",
    "ACGGnFeatureOnly.mf",
    "ACGnFeatureOnly.mf",
    "ACnFeatureOnly.mf",
    "ACTCnFeatureOnly.mf",
    "ACTGnFeatureOnly.mf",
    "ACTnFeatureOnly.mf",
    "ACTTnFeatureOnly.mf",
    "AGATnFeatureOnly.mf",
    "AGCCnFeatureOnly.mf",
    "AGCGnFeatureOnly.mf",
    "AGCnFeatureOnly.mf",
    "AGCTnFeatureOnly.mf",
    "AGGCnFeatureOnly.mf",
    "AGGGnFeatureOnly.mf",
    "AGGnFeatureOnly.mf",
    "AGGTnFeatureOnly.mf",
    "AGnFeatureOnly.mf",
    "AGTCnFeatureOnly.mf",
    "AGTGnFeatureOnly.mf",
    "AGTnFeatureOnly.mf",
    "AGTTnFeatureOnly.mf",
    "AnFeatureOnly.mf",
    "APhasedRepeatsFeatureOnly.mf",
    "ATCCnFeatureOnly.mf",
    "ATCnFeatureOnly.mf",
    "ATCTnFeatureOnly.mf",
    "ATGCnFeatureOnly.mf",
    "ATGGnFeatureOnly.mf",
    "ATGnFeatureOnly.mf",
    "ATGTnFeatureOnly.mf",
    "ATnFeatureOnly.mf",
    "ATTCnFeatureOnly.mf",
    "ATTGnFeatureOnly.mf",
    "ATTnFeatureOnly.mf",
    "ATTTnFeatureOnly.mf",
    "CCCGnFeatureOnly.mf",
    "CCCTnFeatureOnly.mf",
    "CCGGnFeatureOnly.mf",
    "CCGnFeatureOnly.mf",
    "CCGTnFeatureOnly.mf",
    "CCTGnFeatureOnly.mf",
    "CCTnFeatureOnly.mf",
    "CCTTnFeatureOnly.mf",
    "CGCTnFeatureOnly.mf",
    "CGGGnFeatureOnly.mf",
    "CGGnFeatureOnly.mf",
    "CGGTnFeatureOnly.mf",
    "CGnFeatureOnly.mf",
    "CGTnFeatureOnly.mf",
    "CGTTnFeatureOnly.mf",
    "CnFeatureOnly.mf",
    "CTGGnFeatureOnly.mf",
    "CTGnFeatureOnly.mf",
    "CTGTnFeatureOnly.mf",
    "CTnFeatureOnly.mf",
    "CTTGnFeatureOnly.mf",
    "CTTnFeatureOnly.mf",
    "CTTTnFeatureOnly.mf",
    "DirectRepeatsFeatureOnly.mf",
    "EmptyFeatureOnly.mf",
    "GGGTnFeatureOnly.mf",
    "GGTnFeatureOnly.mf",
    "GGTTnFeatureOnly.mf",
    "GnFeatureOnly.mf",
    "GQuadMinusFeatureOnly.mf",
    "GQuadPlusFeatureOnly.mf",
    "GTnFeatureOnly.mf",
    "GTTnFeatureOnly.mf",
    "GTTTnFeatureOnly.mf",
    "InvertedRepeatsFeatureOnly.mf",
    "MirrorRepeatsFeatureOnly.mf",
    "TnFeatureOnly.mf",
    "ZDNAMotifsFeatureOnly.mf"
  )
listEmpty <-
  c(
    "AAACnFeatureOnly.mfEmptyTmp",
    "AAAGnFeatureOnly.mfEmptyTmp",
    "AAATnFeatureOnly.mfEmptyTmp",
    "AACCnFeatureOnly.mfEmptyTmp",
    "AACGnFeatureOnly.mfEmptyTmp",
    "AACnFeatureOnly.mfEmptyTmp",
    "AACTnFeatureOnly.mfEmptyTmp",
    "AAGCnFeatureOnly.mfEmptyTmp",
    "AAGGnFeatureOnly.mfEmptyTmp",
    "AAGnFeatureOnly.mfEmptyTmp",
    "AAGTnFeatureOnly.mfEmptyTmp",
    "AATCnFeatureOnly.mfEmptyTmp",
    "AATGnFeatureOnly.mfEmptyTmp",
    "AATnFeatureOnly.mfEmptyTmp",
    "AATTnFeatureOnly.mfEmptyTmp",
    "ACAGnFeatureOnly.mfEmptyTmp",
    "ACATnFeatureOnly.mfEmptyTmp",
    "ACCCnFeatureOnly.mfEmptyTmp",
    "ACCGnFeatureOnly.mfEmptyTmp",
    "ACCnFeatureOnly.mfEmptyTmp",
    "ACCTnFeatureOnly.mfEmptyTmp",
    "ACGGnFeatureOnly.mfEmptyTmp",
    "ACGnFeatureOnly.mfEmptyTmp",
    "ACnFeatureOnly.mfEmptyTmp",
    "ACTCnFeatureOnly.mfEmptyTmp",
    "ACTGnFeatureOnly.mfEmptyTmp",
    "ACTnFeatureOnly.mfEmptyTmp",
    "ACTTnFeatureOnly.mfEmptyTmp",
    "AGATnFeatureOnly.mfEmptyTmp",
    "AGCCnFeatureOnly.mfEmptyTmp",
    "AGCGnFeatureOnly.mfEmptyTmp",
    "AGCnFeatureOnly.mfEmptyTmp",
    "AGCTnFeatureOnly.mfEmptyTmp",
    "AGGCnFeatureOnly.mfEmptyTmp",
    "AGGGnFeatureOnly.mfEmptyTmp",
    "AGGnFeatureOnly.mfEmptyTmp",
    "AGGTnFeatureOnly.mfEmptyTmp",
    "AGnFeatureOnly.mfEmptyTmp",
    "AGTCnFeatureOnly.mfEmptyTmp",
    "AGTGnFeatureOnly.mfEmptyTmp",
    "AGTnFeatureOnly.mfEmptyTmp",
    "AGTTnFeatureOnly.mfEmptyTmp",
    "AnFeatureOnly.mfEmptyTmp",
    "APhasedRepeatsFeatureOnly.mfEmptyTmp",
    "ATCCnFeatureOnly.mfEmptyTmp",
    "ATCnFeatureOnly.mfEmptyTmp",
    "ATCTnFeatureOnly.mfEmptyTmp",
    "ATGCnFeatureOnly.mfEmptyTmp",
    "ATGGnFeatureOnly.mfEmptyTmp",
    "ATGnFeatureOnly.mfEmptyTmp",
    "ATGTnFeatureOnly.mfEmptyTmp",
    "ATnFeatureOnly.mfEmptyTmp",
    "ATTCnFeatureOnly.mfEmptyTmp",
    "ATTGnFeatureOnly.mfEmptyTmp",
    "ATTnFeatureOnly.mfEmptyTmp",
    "ATTTnFeatureOnly.mfEmptyTmp",
    "CCCGnFeatureOnly.mfEmptyTmp",
    "CCCTnFeatureOnly.mfEmptyTmp",
    "CCGGnFeatureOnly.mfEmptyTmp",
    "CCGnFeatureOnly.mfEmptyTmp",
    "CCGTnFeatureOnly.mfEmptyTmp",
    "CCTGnFeatureOnly.mfEmptyTmp",
    "CCTnFeatureOnly.mfEmptyTmp",
    "CCTTnFeatureOnly.mfEmptyTmp",
    "CGCTnFeatureOnly.mfEmptyTmp",
    "CGGGnFeatureOnly.mfEmptyTmp",
    "CGGnFeatureOnly.mfEmptyTmp",
    "CGGTnFeatureOnly.mfEmptyTmp",
    "CGnFeatureOnly.mfEmptyTmp",
    "CGTnFeatureOnly.mfEmptyTmp",
    "CGTTnFeatureOnly.mfEmptyTmp",
    "CnFeatureOnly.mfEmptyTmp",
    "CTGGnFeatureOnly.mfEmptyTmp",
    "CTGnFeatureOnly.mfEmptyTmp",
    "CTGTnFeatureOnly.mfEmptyTmp",
    "CTnFeatureOnly.mfEmptyTmp",
    "CTTGnFeatureOnly.mfEmptyTmp",
    "CTTnFeatureOnly.mfEmptyTmp",
    "CTTTnFeatureOnly.mfEmptyTmp",
    "DirectRepeatsFeatureOnly.mfEmptyTmp",
    "EmptyFeatureOnly.mfEmptyTmp",
    "GGGTnFeatureOnly.mfEmptyTmp",
    "GGTnFeatureOnly.mfEmptyTmp",
    "GGTTnFeatureOnly.mfEmptyTmp",
    "GnFeatureOnly.mfEmptyTmp",
    "GQuadMinusFeatureOnly.mfEmptyTmp",
    "GQuadPlusFeatureOnly.mfEmptyTmp",
    "GTnFeatureOnly.mfEmptyTmp",
    "GTTnFeatureOnly.mfEmptyTmp",
    "GTTTnFeatureOnly.mfEmptyTmp",
    "InvertedRepeatsFeatureOnly.mfEmptyTmp",
    "MirrorRepeatsFeatureOnly.mfEmptyTmp",
    "TnFeatureOnly.mfEmptyTmp",
    "ZDNAMotifsFeatureOnly.mfEmptyTmp"
  )


list <- c(listMotifs, listEmpty)

library(parallel)

# Calculate the number of cores
no_cores <- 8
# Initiate cluster
cl <- makeCluster(no_cores, type = "FORK")

print("PARALLEL")
ptm <- proc.time()
parLapply(cl, list,
          function(motif)
            processMotif(motif))
stopCluster(cl)
proc.time() - ptm


print("Done.")
traceback()