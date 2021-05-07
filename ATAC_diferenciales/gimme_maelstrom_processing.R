# Transformations applied to peak counts prior to gimme maelstrom
library(limma)

merge_mean <- function(tmp){
  # Compute the mean of each pair of replicates
  # Author: areyes
  
  ts_data <- data.frame(row.names = rownames(tmp))
  for (col in 1:(ncol(tmp)/2)) {
    
    ts_data <- cbind(ts_data, data.frame((tmp[,2*col - 1] + tmp[,2*col])/2))
  }
  colnames(ts_data) <- c("8w", "16w", "24w")
  return(ts_data)
}


for (tissue in c("FB", "HB", "MB")) {

  # Read counts 
  file_counts <- paste0("/home/ska/areyes/nfur/ATACseq/counts/atacpeaks_counts_", tissue, ".txt")
  countdata <- read.table(file_counts, header=FALSE, sep="\t", quote = "", row.names = 1)[,-c(2,3)]
  colnames(countdata) <- c("8w_1", "8w_2", "16w_1", "16w_2", "24w_1", "24w_2")
  
  # Log-transform counts
  log_data <- log(countdata+1)
  
  # Quantile normalisation
  norm_peaks <- normalizeQuantiles(log_data)
  
  # Mean-center per row
  mean_centered <- t(apply(norm_peaks, 1, function(x) x - mean(as.numeric(x))))
  mean_centered <- as.data.frame(mean_centered)
  
  # Save files
  mean_centered$loc <- rownames(mean_centered)
  write.table(mean_centered[, c(7, 1:6)], 
              file = paste0("/home/ska/areyes/nfur/ATACseq/gimme_maelstrom/maelstrom_atacpeaks_counts_", tissue, ".txt"),
              quote = F, sep="\t", row.names = F)
  
  # Filter differential peaks
  peaks_diff <- rbind(read.table(paste0("~/nfur/ATACseq/edger/diff_peaks_down_", tissue, ".bed"))[,1:3],
                      read.table(paste0("~/nfur/ATACseq/edger/diff_peaks_up_", tissue, ".bed"))[,1:3])
  peaks_diff$coor <- paste0(peaks_diff$V1, ":", peaks_diff$V2, "-", peaks_diff$V3)
  
  mean_centered <- mean_centered[mean_centered$loc %in% peaks_diff$coor,]
  write.table(mean_centered[, c(7, 1:6)], 
              file = paste0("/home/ska/areyes/nfur/ATACseq/gimme_maelstrom/maelstrom_diff_atacpeaks_counts_", tissue, ".txt"),
              quote = F, sep="\t", row.names = F)

  # Merge replicates (after normalization across samples)
  merged <- merge_mean(norm_peaks)

  # Mean-center per row
  mean_centered <- t(apply(merged, 1, function(x) x - mean(as.numeric(x))))
  mean_centered <- as.data.frame(mean_centered)
  mean_centered$loc <- rownames(mean_centered)
  mean_centered <- mean_centered[mean_centered$loc %in% peaks_diff$coor,]
  
  # Save files
  write.table(mean_centered[, c(4, 1:3)], 
              file = paste0("/home/ska/areyes/nfur/ATACseq/gimme_maelstrom/maelstrom_diff_atacpeaks_counts_merged_", tissue, ".txt"),
              quote = F, sep="\t", row.names = F)
}

