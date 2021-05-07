#### GeneID to gene name function

nfugene_to_dangene <- function(folder)
  # folder variable must contain files ending "txt" with Nfur gene annotation in 
  # the FIRST COLUMN
{
  library(dplyr)
  gn_all <- read.table("/home/areyes/nfur/genome/Nfur_gene_names.txt")[,c(1,3)]
  gn_unique <- gn_all %>% distinct()
  setwd(folder)
  
  # Remove empty files
  docs <- list.files(pattern = "*.txt") 
  file.remove(docs[file.size(docs) == 0])
  
  for (file in list.files(folder, pattern = "txt")) {
    gid <- read.table(file)[,1]
    gn <- gn_unique[gn_unique[,1] %in% gid,2]
    gn <- gn[!gn %in% c("protein_coding", "pseudo_transcript", "repeat_transcript",
                        "dodgy_transcript")]
    write.table(gn, paste(folder, "/", unlist(strsplit(file, ".txt")), "_genename.txt", sep = ""), 
                quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
}

### GO Terms PANTHER

pantherGO <- function(folder)
{
  library(curl)
  library(jsonlite)
  
  setwd(folder)
  dir.create("GO")
  
  # Remove empty files
  docs <- list.files(pattern = "*genename.txt") 
  file.remove(docs[file.size(docs) == 0])
  
  for (infile in list.files(folder, pattern = "genename")) {
    
    cat(infile, "   processing...\n")
    
    # Only the 1000th first genes of the cluster will be read 
    geneInputList <- paste(as.character(read.table(infile, header = F)[1:1000,1]), collapse = ",")
    
    url <- paste("http://pantherdb.org/services/oai/pantherdb/enrich/overrep?",
                 "geneInputList=", geneInputList, 
                 "&organism=7955", # human=9606
                 "&annotDataSet=GO%3A0008150", # MF=GO:0003674  CC=GO:0005575
                 "&enrichmentTestType=FISHER",
                 "&correction=FDR",
                 sep = "")
    
    outfile <- paste("GO/", head(unlist(strsplit(infile, "genename")), n=1), "GO.txt", sep="")
    
    req <- NULL
    error <- try(req <- curl_fetch_memory(url),silent=T)
    
    if (!is.null(req)){
      content <- prettify(rawToChar(req$content))
      json_content <- fromJSON(content)
      table <- json_content$results$result
      sig_terms <- table[table$fdr<0.05,]
      write.table(sig_terms, outfile, col.names = T, row.names = F, quote = F, sep = "\t")
    }
    else{
      writeLines(error, outfile)
    }
  } 
}


merge_mean <- function(tmp){
  # Compute the mean of each pair of replicates
  
  ts_data <- data.frame(row.names = rownames(tmp))
  for (col in 1:(ncol(tmp)/2)) {
    
    ts_data <- cbind(ts_data, data.frame((tmp[,2*col - 1] + tmp[,2*col])/2))
  }
  return(ts_data)
}

RowVar <- function(x, ...) {
  # Compute the variance of each gene across all time points
  
  rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
}

normalize <- function(x) {
  # Apply MinMax scaling to every feature (column)
  
  return ((x - min(x)) / (max(x) - min(x)))
}


# Similar function to mfuzz.plot2 from Mfuzz package but plot shows number of
# genes beloging to each cluster. Added code at 82 and 125 lines.

mfuzz.plot2.areyes <- function (eset, cl, mfrow = c(1, 1), colo, min.mem = 0, time.labels, 
          time.points, ylim.set = c(0, 0), xlab = "Time", ylab = "Expression changes", 
          x11 = TRUE, ax.col = "black", bg = "white", col.axis = "black", 
          col.lab = "black", col.main = "black", col.sub = "black", 
          col = "black", centre = FALSE, centre.col = "black", centre.lwd = 2, 
          Xwidth = 5, Xheight = 5, single = FALSE, multiple = FALSE, ...) 
{
  clusterindex <- cl[[3]]
  memship <- cl[[4]]
  memship[memship < min.mem] <- -1
  colorindex <- integer(dim(exprs(eset))[[1]])
  genes_list <- acore(eset, cl = cl, min.acore = min.mem)
  n_genes<-unlist(lapply(genes_list, nrow))
  if (missing(colo)) {
    colo <- c("#FF0000", "#FF1800", "#FF3000", "#FF4800", 
              "#FF6000", "#FF7800", "#FF8F00", "#FFA700", "#FFBF00", 
              "#FFD700", "#FFEF00", "#F7FF00", "#DFFF00", "#C7FF00", 
              "#AFFF00", "#97FF00", "#80FF00", "#68FF00", "#50FF00", 
              "#38FF00", "#20FF00", "#08FF00", "#00FF10", "#00FF28", 
              "#00FF40", "#00FF58", "#00FF70", "#00FF87", "#00FF9F", 
              "#00FFB7", "#00FFCF", "#00FFE7", "#00FFFF", "#00E7FF", 
              "#00CFFF", "#00B7FF", "#009FFF", "#0087FF", "#0070FF", 
              "#0058FF", "#0040FF", "#0028FF", "#0010FF", "#0800FF", 
              "#2000FF", "#3800FF", "#5000FF", "#6800FF", "#8000FF", 
              "#9700FF", "#AF00FF", "#C700FF", "#DF00FF", "#F700FF", 
              "#FF00EF", "#FF00D7", "#FF00BF", "#FF00A7", "#FF008F", 
              "#FF0078", "#FF0060", "#FF0048", "#FF0030", "#FF0018")
  }
  else {
    if (colo == "fancy") {
      fancy.blue <- c(c(255:0), rep(0, length(c(255:0))), 
                      rep(0, length(c(255:150))))
      fancy.green <- c(c(0:255), c(255:0), rep(0, length(c(255:150))))
      fancy.red <- c(c(0:255), rep(255, length(c(255:0))), 
                     c(255:150))
      colo <- rgb(b = fancy.blue/255, g = fancy.green/255, 
                  r = fancy.red/255)
    }
  }
  colorseq <- seq(0, 1, length = length(colo))
  if (multiple){
    for (j in multiple) {
      tmp <- exprs(eset)[clusterindex == j, , drop = FALSE]
      tmpmem <- memship[clusterindex == j, j]
      if (((j - 1)%%(mfrow[1] * mfrow[2])) == 0 | single) {
        if (x11) 
          X11(width = Xwidth, height = Xheight)
        if (sum(clusterindex == j) == 0) {
          ymin <- -1
          ymax <- +1
        }
        else {
          ymin <- min(tmp)
          ymax <- max(tmp)
        }
        if (sum(ylim.set == c(0, 0)) == 2) {
          ylim <- c(ymin, ymax)
        }
        else {
          ylim <- ylim.set
        }
        if (!is.na(sum(mfrow))) {
          par(mfrow = mfrow, bg = bg, col.axis = col.axis, 
              col.lab = col.lab, col.main = col.main, col.sub = col.sub, 
              col = col)
        }
        else {
          par(bg = bg, col.axis = col.axis, col.lab = col.lab, 
              col.main = col.main, col.sub = col.sub, col = col)
        }
        xlim.tmp <- c(1, dim(exprs(eset))[[2]])
        if (!(missing(time.points))) 
          xlim.tmp <- c(min(time.points), max(time.points))
        plot.default(x = NA, xlim = xlim.tmp, ylim = ylim, 
                     xlab = xlab, ylab = ylab, main = paste("Cluster", 
                                                            j), axes = FALSE, ...)
        mtext(paste(cl[[2]][j], "genes"), cex = 0.7)
        if (missing(time.labels) && missing(time.points)) {
          axis(1, 1:dim(exprs(eset))[[2]], c(1:dim(exprs(eset))[[2]]), 
               col = ax.col, ...)
          axis(2, col = ax.col, ...)
        }
        if (missing(time.labels) && !(missing(time.points))) {
          axis(1, time.points, 1:length(time.points), time.points, 
               col = ax.col, ...)
          axis(2, col = ax.col, ...)
        }
        if (missing(time.points) & !(missing(time.labels))) {
          axis(1, 1:dim(exprs(eset))[[2]], time.labels, 
               col = ax.col, ...)
          axis(2, col = ax.col, ...)
        }
        if (!(missing(time.points)) & !(missing(time.labels))) {
          axis(1, time.points, time.labels, col = ax.col, 
               ...)
          axis(2, col = ax.col, ...)
        }
      }
      else {
        if (sum(clusterindex == j) == 0) {
          ymin <- -1
          ymax <- +1
        }
        else {
          ymin <- min(tmp)
          ymax <- max(tmp)
        }
        if (sum(ylim.set == c(0, 0)) == 2) {
          ylim <- c(ymin, ymax)
        }
        else {
          ylim <- ylim.set
        }
        xlim.tmp <- c(1, dim(exprs(eset))[[2]])
        if (!(missing(time.points))) 
          xlim.tmp <- c(min(time.points), max(time.points))
        plot.default(x = NA, xlim = xlim.tmp, ylim = ylim, 
                     xlab = xlab, ylab = ylab, main = paste("Cluster", 
                                                            j), axes = FALSE, ...)
        
        mtext(paste(n_genes[j], "genes"), cex = 0.7)
        #mtext(paste(cl[[2]][j], "genes"), cex = 0.7)

                if (missing(time.labels) && missing(time.points)) {
          axis(1, 1:dim(exprs(eset))[[2]], c(1:dim(exprs(eset))[[2]]), 
               col = ax.col, ...)
          axis(2, col = ax.col, ...)
        }
        if (missing(time.labels) && !(missing(time.points))) {
          axis(1, time.points, 1:length(time.points), time.points, 
               col = ax.col, ...)
          axis(2, col = ax.col, ...)
        }
        if (missing(time.points) & !(missing(time.labels))) {
          axis(1, 1:dim(exprs(eset))[[2]], time.labels, 
               col = ax.col, ...)
          axis(2, col = ax.col, ...)
        }
        if (!(missing(time.points)) & !(missing(time.labels))) {
          axis(1, time.points, time.labels, col = ax.col, 
               ...)
          axis(2, col = ax.col, ...)
        }
      }
      if (length(tmpmem) > 0) {
        for (jj in 1:(length(colorseq) - 1)) {
          tmpcol <- (tmpmem >= colorseq[jj] & tmpmem <= 
                       colorseq[jj + 1])
          if (sum(tmpcol) > 0) {
            tmpind <- which(tmpcol)
            for (k in 1:length(tmpind)) {
              if (missing(time.points)) {
                lines(tmp[tmpind[k], ], col = colo[jj])
              }
              else lines(time.points, tmp[tmpind[k], ], 
                         col = colo[jj])
            }
          }
        }
      }
      if (centre) {
        lines(cl[[1]][j, ], col = centre.col, lwd = centre.lwd)
      }
      if (single) 
        return()
    }
  }
  else{
    for (j in 1:dim(cl[[1]])[[1]]) {
      if (single) 
        j <- single
      tmp <- exprs(eset)[clusterindex == j, , drop = FALSE]
      tmpmem <- memship[clusterindex == j, j]
      if (((j - 1)%%(mfrow[1] * mfrow[2])) == 0 | single) {
        if (x11) 
          X11(width = Xwidth, height = Xheight)
        if (sum(clusterindex == j) == 0) {
          ymin <- -1
          ymax <- +1
        }
        else {
          ymin <- min(tmp)
          ymax <- max(tmp)
        }
        if (sum(ylim.set == c(0, 0)) == 2) {
          ylim <- c(ymin, ymax)
        }
        else {
          ylim <- ylim.set
        }
        if (!is.na(sum(mfrow))) {
          par(mfrow = mfrow, bg = bg, col.axis = col.axis, 
              col.lab = col.lab, col.main = col.main, col.sub = col.sub, 
              col = col)
        }
        else {
          par(bg = bg, col.axis = col.axis, col.lab = col.lab, 
              col.main = col.main, col.sub = col.sub, col = col)
        }
        xlim.tmp <- c(1, dim(exprs(eset))[[2]])
        if (!(missing(time.points))) 
          xlim.tmp <- c(min(time.points), max(time.points))
        plot.default(x = NA, xlim = xlim.tmp, ylim = ylim, 
                     xlab = xlab, ylab = ylab, main = paste("Cluster", 
                                                            j), axes = FALSE, ...)
        mtext(paste(cl[[2]][j], "genes"), cex = 0.7)
        if (missing(time.labels) && missing(time.points)) {
          axis(1, 1:dim(exprs(eset))[[2]], c(1:dim(exprs(eset))[[2]]), 
               col = ax.col, ...)
          axis(2, col = ax.col, ...)
        }
        if (missing(time.labels) && !(missing(time.points))) {
          axis(1, time.points, 1:length(time.points), time.points, 
               col = ax.col, ...)
          axis(2, col = ax.col, ...)
        }
        if (missing(time.points) & !(missing(time.labels))) {
          axis(1, 1:dim(exprs(eset))[[2]], time.labels, 
               col = ax.col, ...)
          axis(2, col = ax.col, ...)
        }
        if (!(missing(time.points)) & !(missing(time.labels))) {
          axis(1, time.points, time.labels, col = ax.col, 
               ...)
          axis(2, col = ax.col, ...)
        }
      }
      else {
        if (sum(clusterindex == j) == 0) {
          ymin <- -1
          ymax <- +1
        }
        else {
          ymin <- min(tmp)
          ymax <- max(tmp)
        }
        if (sum(ylim.set == c(0, 0)) == 2) {
          ylim <- c(ymin, ymax)
        }
        else {
          ylim <- ylim.set
        }
        xlim.tmp <- c(1, dim(exprs(eset))[[2]])
        if (!(missing(time.points))) 
          xlim.tmp <- c(min(time.points), max(time.points))
        plot.default(x = NA, xlim = xlim.tmp, ylim = ylim, 
                     xlab = xlab, ylab = ylab, main = paste("Cluster", 
                                                            j), axes = FALSE, ...)
        
        mtext(paste(n_genes[j], "genes"), cex = 0.7)
        #mtext(paste(cl[[2]][j], "genes"), cex = 0.7)
        
        if (missing(time.labels) && missing(time.points)) {
          axis(1, 1:dim(exprs(eset))[[2]], c(1:dim(exprs(eset))[[2]]), 
               col = ax.col, ...)
          axis(2, col = ax.col, ...)
        }
        if (missing(time.labels) && !(missing(time.points))) {
          axis(1, time.points, 1:length(time.points), time.points, 
               col = ax.col, ...)
          axis(2, col = ax.col, ...)
        }
        if (missing(time.points) & !(missing(time.labels))) {
          axis(1, 1:dim(exprs(eset))[[2]], time.labels, 
               col = ax.col, ...)
          axis(2, col = ax.col, ...)
        }
        if (!(missing(time.points)) & !(missing(time.labels))) {
          axis(1, time.points, time.labels, col = ax.col, 
               ...)
          axis(2, col = ax.col, ...)
        }
      }
      if (length(tmpmem) > 0) {
        for (jj in 1:(length(colorseq) - 1)) {
          tmpcol <- (tmpmem >= colorseq[jj] & tmpmem <= 
                       colorseq[jj + 1])
          if (sum(tmpcol) > 0) {
            tmpind <- which(tmpcol)
            for (k in 1:length(tmpind)) {
              if (missing(time.points)) {
                lines(tmp[tmpind[k], ], col = colo[jj])
              }
              else lines(time.points, tmp[tmpind[k], ], 
                         col = colo[jj])
            }
          }
        }
      }
      if (centre) {
        lines(cl[[1]][j, ], col = centre.col, lwd = centre.lwd)
      }
      if (single) 
        return()
    }
  }
  
}
