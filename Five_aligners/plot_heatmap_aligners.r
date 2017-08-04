rescale <- function(x)
  (x) / (max(x)) * 100

getColor <- function(significance, sign) {
  if (significance < 1e-5) {
    significance <- 1e-5
  }
  #significance[significance>0.05]=1 # threshold at 5%
  log_pvalue = -log10(significance)
  log_pvalue = sign * log_pvalue     # sign_difference -1 for negative differences, +1 for   positive differences
  colfunc = colorRampPalette(c("navy", "white", "red"))(n = 100 - 1)
  # 99 colors from blue to white to red.
  # Blue will correspond to log_pvalue=-5
  # white log_pvalue=0
  # red log_pvalue=5
  color_limits = seq (-5, 5, length.out = 100)
  col = c(colfunc[unlist(lapply(log_pvalue, function(log)
    which((log >= color_limits[1:99]) &
            (log <= color_limits[2:100])
    )[1]))])
  if (significance > 0.05) {
    col <- "white" #non-significant not colored
  }
  return(col)
}

getParametersForTechnology <- function(technology) {
  print(technology[, 1])
  decimal_points <- 4
  columns<-8
  signs <- as.numeric(as.matrix(technology[seq(6, (7*columns)+1, columns)]))
  significance <- as.numeric(as.matrix(technology[seq(7, (7*columns)+1, columns)]))
  effect <-
    round(as.numeric(as.matrix(technology[seq(8, (7*columns)+1, columns)])), decimal_points)
  zeros <-
    round(as.numeric(as.matrix(technology[seq(9, (7*columns)+1, columns)])), decimal_points)
  
  mean_features <- as.numeric(as.matrix(technology[seq(2, (7*columns)+1, columns)]))
  mean_controls <- as.numeric(as.matrix(technology[seq(3, (7*columns)+1, columns)]))
  
  real_perc_diff <- round(as.numeric(mean_features-mean_controls)/mean_controls, decimal_points)

  #signs[signs==0]<-sign(zeros)[signs==0]
  #todo: how to pick color when the test is based on the proportion of 0s?
  print(paste(significance, effect, zeros))
  return(list(
    effect = effect,
    col = col,
    zeros = zeros,
    significance = significance,
    signs = signs,
    real_perc_diff = real_perc_diff
  ))
}

options(scipen=10)
#PLOT HEATMAP

data <-
  as.data.frame(
    read.table(
      "aligners.raw_data_for_Figure_4.txt",
      sep = "\t"
    ),
    optional = TRUE
  )

dataSections <- 3  # MISMATCHES INSERTIONS DELETIONS
rowsPerSection <- (nrow(data)-(dataSections-1)) / dataSections
numAligners <- rowsPerSection - 3
alignerNames <- as.character(data[4:(4+numAligners-1),1])

reorderMojo = c(1,2,5,6,7,3,4)
motifNames <- unique(as.character(unlist(data[2,])))
motifNames <- motifNames[2:length(motifNames)]
motifNames <- gsub(".features.plus.collapsed", "", motifNames)


for (dataSection in 1:dataSections) {
  sectionStart <- 1 + (dataSection-1) * (rowsPerSection+1)
  errorTypeName <- as.character(data[sectionStart,1])

  quartz()
  par(bg = "white", las = 1, pty = "s")
  par(xpd=TRUE)
  
  barplot(
    c(0),
    c(0),
    xaxs = "i",
    yaxs = "i",
    col = colours()[1:6],
    ylim = c(0, 700),
    xlim = c(0, 700),
    bg = "white",
    axes = FALSE
  )
  
  text(seq(50,700,100),800,labels=motifNames[reorderMojo],srt=90,adj=.4)
  
  mtext(c(paste0(alignerNames," "),rep("",7)),
        side = 2,
        at = seq(640, 0, -100))

  sectionDataStart <- sectionStart+3
  for (alignerNum in (1:numAligners)) {
    technology <- alignerNum + (sectionDataStart-1)
    print(technology)
    result = getParametersForTechnology(data[technology,])
    print(result)
    
    i <- 0
    for (x in seq(0, 699, 100)) {
      i <- i + 1
      iFeature <- reorderMojo[i]
  
      y <- 700 - (alignerNum * 100) #skip header
      size <- rescale(abs(result$effect))[iFeature]
      zero <- rescale(abs(result$zeros))[iFeature]
      sizeToPlot <- 100
      
      backgroundColor1 <- "white"
      backgroundColor2 <- "white"
      
      if (!is.na(result$significance[iFeature])) {
        #plot just rectangle
        polygon(c(x, x,                x + (sizeToPlot), x + (sizeToPlot)),
                c(y, y + (sizeToPlot), y + (sizeToPlot), y),
                col = getColor(result$significance[iFeature], result$signs[iFeature]), lty=0)
      }
        
      cell_value<-result$real_perc_diff[iFeature]*100
      if (!is.na(cell_value)) {
        cell_value<-round(cell_value)
        if (cell_value==0) {
          cell_value<-paste0("~",cell_value,"%")
        } else {
          cell_value<-paste0(cell_value,"%")
        }
      }
      text(x + 50,
           y + 50,
           labels = cell_value,
           col = ifelse(result$significance[iFeature] >0.05, "black","white"),
           cex=1.25
        )
    }
  }
  grid(col = "black", lty = 1)
  text(0,700-(numAligners*100)-20, errorTypeName)
}
