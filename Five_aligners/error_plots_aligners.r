options(scipen = 10)
options(error = traceback)
options(show.error.locations = TRUE)
#par(
#  mfrow = c(4, 4),
#  cex.main = 1.3,
#  cex.sub = 1.4,
#  cex.axis = 1.2
#)
require("plyr")
require("dplyr")
require("hash")
require("data.table")

testCoordinates <- function(data, dataEmpty) {
  feature_lengths <- data[, 3] - data[, 2]
  control_lengths <- dataEmpty[, 3] - dataEmpty[, 2]
  return(all.equal(feature_lengths, control_lengths))
}

getPvalue <- function(grp.1, grp.2, file_id, i, ITER, paired, parametric=FALSE) {
  failure_result <- list(
      pvalue = NA,
      effectSize = NA,
      effectSize_perc = NA,
      nonZeros_perc = NA,
      param_pval = NA
    )

  if ((length(which(!is.na(grp.1))) == 0) || (length(which(!is.na(grp.2))) == 0)) {
    print('group 1 or 2 is \'empty\'')
    return(failure_result)
  }

  print(paste("mean(grp.1) =",mean(grp.1)))
  print(paste("mean(grp.2) =",mean(grp.2)))

  #populate hashes
  hash_name <- i
  tmp_data <- get(hash_name)
  tmp_data[file_id] = grp.1 #features
  assign(hash_name, tmp_data)
  
  #populate hash with matching controls
  hash_name <- paste0(i, "controls")
  tmp_dataEmpty <- get(hash_name)
  tmp_dataEmpty[file_id] = grp.2 #controls
  assign(hash_name, tmp_dataEmpty)
  
  feature_non_zero <-
    round((length(grp.1) - (table(grp.1)["0"])) / length(grp.1), digits = 4)
  control_non_zero <-
    round((length(grp.2) - (table(grp.2)["0"])) / length(grp.2), digits = 4)
  
  if (is.na(feature_non_zero)) {
    feature_non_zero <- 1 #no zeros exist, all values are non-zero
  }
  
  if (is.na(control_non_zero)) {
    control_non_zero <- 1 #no zeros exist, all values are non-zero
  }
  
  print(
    paste(
      "feature_non_zero: ",
      feature_non_zero,
      "; control_non_zero:",
      control_non_zero
    )
  )
  nonZeros_perc <-
    as.numeric((feature_non_zero - control_non_zero) / control_non_zero)
  
  # Function for calculating two-part statistics
  # Modified from Taylor, Pollard (2009) Hypothesis Tests for Point-Mass Mixture
  # Data with Application to `Omics Data with Many Zero Values
  TwoPart <- function(data,
                      group,
                      test = "t.test",
                      point.mass = 0) {
    Index1 <- c(group == 1)
    Group1 <- data[Index1]
    Group0 <- data[!Index1]
    n1 <- length(Group1)
    n2 <- length(Group0)
    obs <- c(n1, n2)
    success <- c(sum(Group1 != point.mass), sum(Group0 != point.mass))
    pointmass <- obs - success
    if (sum(success) == 0) {
      T2 <- 0
      difference = 0
      difference_perc = 0
      B2 <- 0
    } else if ((success[1] == 0) | (success[2] == 0)) {
      T2 <- 0
      difference = NA
      difference_perc = NA
      B2 <- prop.test(pointmass, obs)$statistic
    } else if ((success[1] == 1) | (success[2] == 1)) {
      T2 <- 0
      difference = NA
      difference_perc = NA
      B2 <- prop.test(pointmass, obs)$statistic
    } else {
      uniq1 <- length(unique(Group1[Group1 != point.mass]))
      uniq2 <- length(unique(Group0[Group0 != point.mass]))
      if ((uniq1 < 2) & (uniq2 < 2)) {
        T2 <- 0
        difference = NA
        difference_perc = NA
        if (sum(pointmass) == 0) {
          B2 <- 0
        } else {
          B2 <- prop.test(pointmass, obs)$statistic
        }
      } else if (sum(pointmass) == 0) {
        B2 <- 0
        if (test == "t.test") {
          T2 <- Reduce(c, by(data, group, mean, na.rm = TRUE))
          difference <- T2[2] - T2[1]
          difference_perc <- difference / T2[1]
          T2 <- t.test(data~group)$statistic^2
        }
      } else {
        B2 <- prop.test(pointmass, obs)$statistic
        contIndex <- data != point.mass
        cont <- data[contIndex]
        cGroup <- group[contIndex]
        n1c <- sum(cGroup == 1)
        n2c <- sum(cGroup == 0)
        if (test == "t.test") {
          T2 <- Reduce(c, by(cont, cGroup, mean, na.rm = TRUE))
          difference <- T2[2] - T2[1]
          difference_perc <- difference / T2[1]
          T2 <- t.test(cont~cGroup)$statistic^2
        }
      }
    }
    X2 <- B2 + T2
    statistic = X2
    if ((T2==0)|(B2==0)) {
      X2pv <- 1-pchisq(X2,1)
    } else {
      X2pv <- 1-pchisq(X2,2)
    }
    return(
      list(
        statistic = statistic,
        difference = difference,
        difference_perc = difference_perc,
        nonZeros_perc = nonZeros_perc,
        param_pval = X2pv
      )
    )
  }
  
  groups <- c(rep(1, length(grp.1)), rep(0, length(grp.2)))
  data <- c(grp.1, grp.2)
  
  if ((paired) & (length(grp.1) != length(grp.2)))
    stop('Different sample size in the two groups, in a paired test!')
  
  
  #diff(by(data, groups, getMQuantile))
  
  #s<- sample(groups, length(groups), FALSE)
  #diff(by(data, s, getMQuantile))
  
  if(!parametric){
    max.iter <- ITER - 1 #500
    examples <- unlist(lapply(1:max.iter, function(x) {
      if (paired) {
        perm = which(sample(c(TRUE, FALSE), length(grp.1), replace = TRUE))
        groups_perm = groups
        groups_perm[perm] = groups[length(grp.1) + perm]
        groups_perm[length(grp.1) + perm] = groups[perm]
      } else{
        groups_perm = sample(groups)
      }
      #temp=Reduce(cbind,by(data, groups_perm, quantile, probs=c(.05, .25, .5, .75, .95),na.rm=TRUE))
      #temp=sum((temp[,2]-temp[,1])^2)
      temp = TwoPart(data, groups_perm, test = "t.test")$statistic  
      return(temp)
    }))
  }
  
  #test.diff <- Reduce(cbind,by(data, groups, quantile, probs=c(.05, .25, .5, .75, .95),na.rm=TRUE))
  #test.diff=sum((test.diff[,2]-test.diff[,1])^2)
  two.part = TwoPart(data, groups, test = "t.test")   
  test.diff = two.part$statistic
  sign_difference = sign(two.part$difference)
  effectSize = two.part$difference
  effectSize_perc = two.part$difference_perc
  if(parametric){
    pvalue=two.part$param_pval
  }else{
    # two-tailed test
    pvalue <-
      (sum(examples >= test.diff) + 1) / (max.iter + 1) #from http://danielnee.com/tag/p-value/
  }
  
  # RSH: this won't work, because I no longer pass in file
  #hist(
  # examples,
  #  col = ifelse(pvalue < 0.05, 'red', 'blue'),
  #  breaks = 50,
  #  main = paste("RPerm", main = gsub(".txt", "", gsub(
  #    "ERRORS", "", gsub("10000.", "", gsub(".mf", "", file))
  #  ))),
  #  xlab = "",
  #  sub = paste(
  #    "pvalue: ",
  #    round(pvalue, 4),
  #    " | C:",
  #    length(grp.1),
  #    " | T:",
  #    length(grp.2),
  #    sep = ""
  #  )
  #)
  #abline(v = test.diff, col = "black", lwd = 4)
  return(
    list(
      pvalue = pvalue,
      effectSize = effectSize,
      effectSize_perc = effectSize_perc,
      nonZeros_perc = nonZeros_perc
    )
  )
}

error_plots.run <- function(directory, filenames, paired, ITER, parametric=FALSE) {
  print(paste0("reading from ",directory))
  print(paste0("  ",filenames," paired=",paired," ITER=",ITER," parametric=",parametric))
  
  #PLOT
  for (i in c("TOTAL", "MISMATCHES", "INSERTION", "DELETION")) {
    #c("TOTAL", "MISMATCHES", "INSERTION", "DELETION")
    for (file in filenames) {
      print("")
      print(paste("==========",basename(directory),file,i,"=========="))
      file_id <- gsub(".features.collapsed", "", file)   # this is our key for the hashes

      FEATURENAMES[file_id] = T

      data <- read.table(
        file,
        header = F,
        col.names = c(
          "chr",
          "start",
          "end",
          "TOTAL",
          "MISMATCHES",
          "INSERTION",
          "DELETION"
        ),
        fill = TRUE
      )
      
      values <- as.numeric(as.vector(data[[i]]))
      #values <- na.omit(values)
      #print(mean(values, na.rm = TRUE))
      
      #hist(
      #  values,
      #  col = "gold",
      #  xlab = "% error",
      #  main = gsub(".txt", "", gsub(
      #    "ERRORS", "", gsub("10000.", "", gsub(".mf", "", file))
      #  )),
      #  xlim = c(0, 5),
      #  sub = paste(
      #    "n=",
      #    length(values),
      #    " errorRate=",
      #    round(mean(values), 4),
      #    sep =
      #      ""
      #  )
      #)
      
      #load empty windows
      emptyFile <- gsub(".features", ".controls", file)
      emptyFile <- gsub(".chr", ".for_chr", emptyFile)

      if (file.exists(emptyFile)) {
        print(paste("empty file exists, will use matching file", emptyFile))
        
        #individial empty files
        
        dataEmpty <-
          read.table(
            emptyFile,
            header = F,
            col.names = c(
              "chr",
              "start",
              "end",
              "TOTAL",
              "MISMATCHES",
              "INSERTION",
              "DELETION"
            ),
            fill = TRUE
          )
        
        valuesEmpty <-
          as.numeric(as.character(as.vector(dataEmpty[[i]])))
        #valuesEmpty <-na.omit(valuesEmpty)
        #print(mean(valuesEmpty, na.rm = TRUE))
        
        
      } else {
        # RSH: I don't have Controls.mf
        missingFile <- emptyFile
        emptyFile <- "Controls.mf"
        print(paste(
          "empty file,missingFile,DOES NOT exist, will use DEFAULT control file: ",
          emptyFile
        ))
        dataEmpty <-
          read.table(
            emptyFile,
            header = F,
            col.names = c(
              "chr",
              "start",
              "end",
              "TOTAL",
              "MISMATCHES",
              "INSERTION",
              "DELETION"
            ),
            fill = TRUE
          )
        
        valuesEmpty <-
          as.numeric(as.character(as.vector(dataEmpty[[i]])))
        #valuesEmpty <-na.omit(valuesEmpty)
        print(mean(valuesEmpty, na.rm = TRUE))
        
        #populate hash with matching controls
        hash_name <- paste0(i, "controls")
        tmp_dataEmpty <- get(hash_name)
        tmp_dataEmpty[file_id] = valuesEmpty
        assign(hash_name, tmp_dataEmpty)
      }
      
      print(paste("paired status: ", paired))
      
      if (paired == TRUE) {
        #PAIRED DATA
        df <- as.data.frame(cbind(valuesEmpty, values))
        missing <- which(is.na(df$values) |
                           is.na(df$valuesEmpty))
        filled <- !(1:(length(values)) %in% missing)
        result = getPvalue(values[filled],
                           valuesEmpty[filled],
                           file_id,
                           i,
                           ITER,
                           paired, #should be TRUE for paired test when using features of variable length
                           parametric)
        print(result)
        print(
          paste(
            "control sample size without NA: ",
            length(na.omit(valuesEmpty[filled])),
            "feature sample size without NA:",
            length(na.omit(values[filled]))
          )
        )
        print("validate that controls and features have matching lengths")
        #VALIDATE THAT CONTROLS AND FEATURES HAVE MATCHING LENGTHS
        comparison <-
          testCoordinates(data[filled,], dataEmpty[filled,])
        if (comparison == TRUE) {
          print("Features and controls have equal feature lengths. GOOD.")
        } else {
          print("Features and controls have unequal feature lengths. ERROR.")
          warning("Features and controls have unequal feature lengths. ERROR.")
        }
      } else {
        #SINGLE CONTROL
        result = getPvalue(na.omit(values),
                           na.omit(valuesEmpty),
                           file_id,
                           i,
                           ITER,
                           paired) #should be FALSE when using single control
        print(
          paste(
            "control sample size without NA: ",
            length(na.omit(valuesEmpty)),
            "feature sample size without NA:",
            length(na.omit(values))
          )
        )
      }
      
      pv <- as.numeric(result$pvalue)
      effect_size <-
        result$effectSize_perc #this uses percentages as opposed to using effectSizeonly
      nonZeros_perc <- result$nonZeros_perc
      
      hash_name <- paste0(i, "pvalue")
      tmp_data <- get(hash_name)
      tmp_data[file_id] = round(pv, 4) #save pvalue and plot
      assign(hash_name, tmp_data)
      
      hash_name <- paste0(i, "effectSize")
      tmp_data <- get(hash_name)
      tmp_data[file_id] = effect_size #save effect size
      assign(hash_name, tmp_data)
      
      #populate hash with number of non_zeros
      hash_name <- paste0(i, "nonzeros")
      tmp_data <- get(hash_name)
      tmp_data[file_id] = nonZeros_perc #difference in non-zero valus in percentages
      assign(hash_name, tmp_data)
    }
  }
  
  #save.image(file = paste0(basename(getwd()),".renv"))
  
  
  #OUTPUT TABLES
  print("========== PRINTING TABLES ==========")
  names <- names(FEATURENAMES)
  
  if (length(names)<2) {
    # e.g. as.numeric(values(get(paste0(i, "pvalue")))[names]) will fail if names has only one element
    print("The number of motifs to be analyzed is too low. Please add more motifs. ERROR.")
    warning("The number of motifs to be analyzed is too low. Please add more motifs. ERROR.")
  }
  
  #par(mfrow = c(1, 1), cex = 1.4)
  directory <- paste0(getwd(), "/")
  variables <- basename(getwd())
  
  for (i in c("TOTAL", "MISMATCHES", "INSERTION", "DELETION")) {
    print(paste("  ===",basename(directory),i,"==="))
    
    significance <-
      as.numeric(values(get(paste0(i, "pvalue")))[names])
    effectSize <-
      as.numeric(values(get(paste0(i, "effectSize")))[names])
    
    featuresVector <- values(get(paste0(i)))[names]
    meanFeatures <-
      lapply(featuresVector, mean)
    meanNon0Features <-
      vapply(featuresVector, function(x)
        mean(x[x != 0]), numeric(1))
    
    controlsVector <- values(get(paste0(i, "controls")))[names]
    meanControls <-
      lapply(controlsVector, mean)
    meanNon0Controls <-
      vapply(controlsVector, function(x)
        mean(x[x != 0]), numeric(1))
    
    print(
      paste(
        " ",
        names,
        "meanF:",
        meanFeatures,
        "meanC:",
        meanControls,
        "meanNon0F:",
        meanNon0Features,
        "meanNon0C:",
        meanNon0Controls
      )
    )
    
    signs <- sign(as.numeric(meanFeatures) - as.numeric(meanControls))
    nonzeros <- values(get(paste0(i, "nonzeros")))[names]
    
    results_headers <-
      c(
        "mean_features",
        "mean_controls",
        "nonZeros_mean_features",
        "nonZeros_mean_controls",
        "sign",
        "p-value",
        "effect-size",
        "zero-effect"
      )
    
    results <-
      paste(
        meanFeatures,
        meanControls,
        meanNon0Features,
        meanNon0Controls,
        signs,
        significance,
        effectSize,
        nonzeros,
        sep = "\t"
      )

    results_one_line <- paste(results,collapse="\t")
    fn <- paste0(directory, variables, "_", i, ".txt")
    print(paste(" ",fn))
    
    if (file.exists(fn)) {
      print(paste("File",fn,"already exists and will be removed before writing to it."))
      file.remove(fn) #remove file if exists
    }

    write.table(paste(rep(names,each=8),collapse="\t"),
                fn, append = FALSE,
                sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    headers <- paste(results_headers,collapse="\t")
    #write headers
    write.table(paste(rep(headers,length(names)),collapse="\t"),
                fn, append = TRUE,
                sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    #append actual results
    write.table(results_one_line,
                fn, append = TRUE,
                sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
  
}

set.seed(3)

# the current working directory must have these as subdirectories:
aligners <- c("novoalign.extracted","novoalign","bwa","bowtie","last","stampy")

directories <- unlist(lapply(aligners, function(aligner) {
  basedir <- getwd()
  return(paste0(basedir,"/",aligner))
}))

basedir <- getwd()
for (directory in directories) {
  setwd(directory)

  FEATURENAMES <- hash()
  
  TOTAL <- hash()
  INSERTION <- hash()
  DELETION <- hash()
  MISMATCHES <- hash()
  
  TOTALcontrols <- hash()
  INSERTIONcontrols <- hash()
  DELETIONcontrols <- hash()
  MISMATCHEScontrols <- hash()
  
  TOTALpvalue <- hash()
  INSERTIONpvalue <- hash()
  DELETIONpvalue <- hash()
  MISMATCHESpvalue <- hash()
  
  TOTALnonzeros <- hash()
  INSERTIONnonzeros <- hash()
  DELETIONnonzeros <- hash()
  MISMATCHESnonzeros <- hash()
  
  TOTALeffectSize <- hash()
  INSERTIONeffectSize <- hash()
  DELETIONeffectSize <- hash()
  MISMATCHESeffectSize <- hash()
  
  filenames <- list.files(pattern = "*.collapsed", full.names = FALSE)
  filenames <- filenames[grep("controls", filenames, invert = TRUE)]

  paired <- TRUE
  ITER <- 1000
  parametric=TRUE # if FALSE, permutational pvalue, otherwise approximate parametric pvalue
  
  error_plots.run(directory, filenames, paired, ITER, parametric)
  #traceback()

}
setwd(basedir)