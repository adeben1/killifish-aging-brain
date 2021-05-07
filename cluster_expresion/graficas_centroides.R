root.name <- "/home/areyes/nfur/RNAseq/mfuzz results/sophisticated/go_metric/differential_only/9_9_90.3mem/"

# Read centroid files
fb_centroids <- read.table(paste0(root.name, "fb9cl_0.3mem/centroids.txt"))
mb_centroids <- read.table(paste0(root.name, "mb9cl_0.3mem/centroids.txt"))
hb_centroids <- read.table(paste0(root.name, "hb9cl_0.3mem/centroids.txt"))

# Define expression profiles (based on heatmap_corr.pdf)
# write in which expression profile the 1th,2nd... cluster is
# vector digits content can be from 1 to number of expression profiles
fb <- c(3, 9, 7, 5, 7, 4, 1, 8, 3)
mb <- c(4, 9, 8, 5, 7, 1, 6, 10, 3)
hb <- c(9, 8, 7, 5, 3, 2, 6, 10, 1)

time <- c(8, 12, 16, 20, 24)

# If lines() doesn't find any vector it prints nothing

for (i in 1:10) {
  png(file=paste0("/home/areyes/nfur/RNAseq/mfuzz results/sophisticated/go_metric/differential_only/9_9_90.3mem/expression_profiles/", i, ".png"),
      width=500, height=450)
  
  par(mar=c(6.1, 7.1, 4.1, 2.1), mgp=c(3, 1.5, 0), las=1)
  plot(time, fb_centroids[fb==i,2:6][1,], pch=1, lwd=7, col="hotpink", xlab="", 
       ylab="", type="b", xaxt="n", yaxt="n", ylim=c(-1.5, 1.5))
  lines(time, fb_centroids[fb==i, 2:6][2,], xlab="", ylab="", type="b", pch=1, lwd=7, col="hotpink")
  lines(time, mb_centroids[mb==i, 2:6][1,], pch=1, lwd=7, col="forestgreen", xlab="", ylab="", type="b")
  lines(time, mb_centroids[mb==i, 2:6][2,], pch=1, lwd=7, col="forestgreen", xlab="", ylab="", type="b")
  lines(time, hb_centroids[hb==i, 2:6][1,], pch=1, lwd=7, col="purple", xlab="", ylab="", type="b")
  lines(time, hb_centroids[hb==i, 2:6][2,], pch=1, lwd=7, col="purple", xlab="", ylab="", type="b")
  axis(1, at=time, labels = time, cex.axis=2, lwd=0, lwd.ticks = 1)
  axis(2, cex.axis=2, at=c(-1.5, 0, 1.5), )
  title(ylab="Standardized TPMs", line=5, cex.lab=2.5)
  title(xlab="Age (weeks)", line=4, cex.lab=2.5)
  title(main=paste(i), line=1, cex.main=5, font.main=2)
  dev.off()
}

dev.off()
par(cex=4)
par(mar=c(1,1,1,1))
plot(1:2, 1:2, xlab = "", ylab = "", axes = F)
legend("center", c("Fore brain", "Mid brain", "Hind brain"), fill = c("hotpink", "forestgreen", "purple"))

