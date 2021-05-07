##### PUll libraries
library(Mfuzz)
library(Biobase)
library(limma)
library(yarrr)
library(ggplot2)
library(ggh4x)
library(pheatmap)
library(dendextend)
source("/home/areyes/nfur/RNAseq/scripts/aux_functions.R")

cl_number <- c(9, 9, 9) # F M H
size_cluster_plots <- c(3,5) # rows cols
mem <- 0.3 # membership filter

# Select patterns to extract from GO term labels and give a nice description
term_pattern <- c("transport.*ion|ion.*transport", "develop","diff", "neuro", "synap", "cell cycle|divis|spindle", "DNA", "RNA", "ribo")
term_description <- c("Ion\ntrans.", "Dev.","Diff.","Neu.", "Syn.", "Cell\ndiv.", "DNA", "RNA", "Ribo.")

folder <- "/home/areyes/nfur/RNAseq/mfuzz results/sophisticated/go_metric/differential_only/"
folder <- paste0(folder, paste(cl_number, collapse = "_"), mem, "mem/")
dir.create(folder)

file_counts <- c("/home/areyes/nfur/RNAseq/counts/rna_ts_counts_fb.tab",
                 "/home/areyes/nfur/RNAseq/counts/rna_ts_counts_mb.tab",
                 "/home/areyes/nfur/RNAseq/counts/rna_ts_counts_hb.tab")

# Significant genes file
file_sig_genes <- "/home/areyes/nfur/RNAseq/deseq results/lrt/sig_genes_alltissues.txt"
sig_genes <- read.table(file_sig_genes)[,1]

# Part 1 - c-means clustering per tissue

fuzzy_clustering <- function(file_counts, c, size_cluster_plots){
  # Read gene counts
  count <- as.matrix(read.table(file_counts[1], header=TRUE, sep="\t", as.is=TRUE))
  tissue <- tail(unlist(strsplit(gsub(".tab", "", file_counts), "_")), 1)
  
  # Filter based on significance 
  filtered <- count[rownames(count) %in% sig_genes,]
  
  # Read gene length
  gene_length <- read.table("/home/areyes/nfur/RNAseq/Nfu_genelengths.tab", header = T)[,c(1,2)]
  gene_length <- gene_length[order(gene_length$gene),]
  gene_length <- gene_length[gene_length[,1] %in% rownames(filtered),]
  diff <- rownames(filtered)[!rownames(filtered) %in% gene_length[,1]]
  filtered <- filtered[!rownames(filtered) %in% diff,]
  
  #### Pre-processing
  
  # TPM (based on https://support.bioconductor.org/p/91218/)
  x <- filtered / gene_length[,2]
  tpm <- t( t(x) * 1e6 / colSums(x) )
  
  # log2( TPM + 2 )
  tpm_log <- ( log(tpm + 1, 2) )
  
  # Sample aggregation by mean
  eset <- merge_mean(tpm_log)
  colnames(eset) <- c("week8", "week12", "week16", "week20", "week24")
  
  # Remove genes with low TPM (when using significant genes do not filter)
  #eset <- eset[rowMeans(eset) > 1,]

  # Creating eset object
  minimalSet <- ExpressionSet(assayData=as.matrix(eset))
  
  # Standardization of data (mean=0 and sd=1)
  eset <- standardise(minimalSet)
  
  # Quantile normalisation
  eset <- normalizeQuantiles(exprs(eset))
  eset <- ExpressionSet(assayData = eset)
  
  # Removing NAs
  eset <- filter.NA(eset, thres=0.25)
  ts_data <- fill.NA(eset,mode="knnw")
  
  # Check data is OK
  pdf(paste0(folder, tissue, "_expression_levels.pdf"))
  data <- data.frame(expression = as.vector(exprs(ts_data)),
                     sample = rep(colnames(ts_data), each=nrow(ts_data)))
  data$sample <- factor(data$sample, levels = colnames(eset))
  pirateplot(formula = expression ~ sample, data = data, point.o = 0, 
             ylab = "", xlab = "Time (weeks)", xaxt = "")
  text(x = c(1:5)+.07, y = par("usr")[3]-.2, labels = c(8, 12, 16, 20, 24), xpd = NA, adj = 1, cex = 1)
  text(x=5.3, y = 2, labels = paste0("n=", nrow(exprs(ts_data))))
  dev.off()
  
  m <- mestimate(ts_data)
  
  set.seed(1)
  cl <- mfuzz(ts_data,c=c,m=m)
  
  # Plot clusters
  plot_name <- paste0(folder, tissue, "_", c, "cl.pdf")
  mfuzz.plot2.areyes(ts_data, cl = cl, mfrow = size_cluster_plots, xlab = "Time (weeks)",
                     time.labels = c(8, 12, 16, 20, 24), Xwidth = 20, Xheight = 15, min.mem = mem)
  dev.copy2pdf(file = plot_name)
  dev.off()
  
  #### Imprimir cada gen con su cluster
  genes_list <- acore(ts_data, cl = cl, min.acore = mem)
  
  dir_name <- paste(folder, tissue, c, "cl_", mem, "mem", sep = "")
  dir.create(dir_name)
  setwd(dir_name)
  
  # Select genes 
  for (cluster in seq(c)){
    file_name <- paste(cluster, ".txt", sep = "")
    select <- order(genes_list[[cluster]]$MEM.SHIP, decreasing = T)
    write.table(genes_list[[cluster]][select,], file_name, quote = F, row.names = F, col.names = F)
  }
  
  # Save centroids and membership
  write.table(cl$centers, paste0(dir_name, "/centroids.txt"), row.names = T, col.names = F, quote = F)
  write.table(cl$membership, paste0(dir_name, "/membership.txt"), row.names = T, col.names = T, quote = F)
  
  # Convert nfu ids to danre ids
  nfugene_to_dangene(dir_name)
  
  # Run GO Analysis
  pantherGO(dir_name)
}

fuzzy_clustering(file_counts[1], cl_number[1], size_cluster_plots)
fuzzy_clustering(file_counts[2], cl_number[2], size_cluster_plots)
fuzzy_clustering(file_counts[3], cl_number[3], size_cluster_plots)

# Save centroids
fb_centroids_file <- paste0(folder, "fb", cl_number[1], "cl_", mem, "mem/", "centroids.txt")
mb_centroids_file <- paste0(folder, "mb", cl_number[2], "cl_", mem, "mem/", "centroids.txt")
hb_centroids_file <- paste0(folder, "hb", cl_number[3], "cl_", mem, "mem/", "centroids.txt")

fb_centroids <- read.table(fb_centroids_file, header = F)
mb_centroids <- read.table(mb_centroids_file, header = F)
hb_centroids <- read.table(hb_centroids_file, header = F)
all <- rbind(fb_centroids, rbind(mb_centroids, hb_centroids))

# Compute correlation of all centroids (even from different tissues)
c <- cor(t(all[2:6]), method = "pearson")
rownames(c) <- c(paste0(1:cl_number[1], "fb"), paste0(1:cl_number[2], "mb"), paste0(1:cl_number[3], "hb"))
colnames(c) <- c(paste0(1:cl_number[1], "fb"), paste0(1:cl_number[2], "mb"), paste0(1:cl_number[3], "hb"))

# Plot heatmap and generate expresion profiles using hierarchical clustering
pheatmap(c, border_color = NA, show_rownames = F, show_colnames = F,fontsize = 15,
         legend = T)

happy = "n"
while (happy == "n") {
  k <- as.integer(readline(prompt = "Enter k divisions: "))
  pheatmap(c, border_color = NA, show_rownames = F, show_colnames = T,fontsize = 15,
           legend = T, cutree_rows = k, cutree_cols = k)
  happy <- readline(prompt = "Are you happy? [y/n] ")
}

d <- dist(c, method = "euclidean")
hclust <- set(as.dendrogram(hclust(d, method = "complete")), "branches_lwd", 2)
a <- cutree_1k.dendrogram(hclust, k)

# Part 3 - Plot similarities among GO terms

# Read GO terms from files and saves term label, cluster and Benjamini
read_GOs <- function(tissue, clusters, folder){
  setwd(folder)
  
  yy <- data.frame()
  
  for (cluster in clusters) {
    
    # Check if cluster n exists (it wont in case there is no gene under mem filter)
    if(file.exists(paste(cluster, "_GO.txt", sep=""))){

      xx <-read.csv(paste(cluster, "_GO.txt", sep=""), header = T, sep = "\t")
      
      if (nrow(xx) == 0) {
        xx[1,] <- c(NA, NA, 1, NA, NA, NA, NA, NA)
        colnames(xx)[7] <- c("term.id")
        xx$plus_minus <- NA
        xx$term.label <- NA
      }
      
      xx$Condition <- as.character(paste0(cluster, tissue))
      xx$fdr <- -log(xx$fdr, 10)
      yy <- rbind(yy, xx)
    }
  }
  
  data.frame("GOs" = yy$term.label, 
             "Condition" = yy$Condition,
             "Benjamini" = yy$fdr)
}

clusters_fb <- 1:cl_number[1]
clusters_hb <- 1:cl_number[3]
clusters_mb <- 1:cl_number[2]

go_fb <- read_GOs("fb", clusters_fb, paste0(folder, "fb", cl_number[1], "cl_", mem, "mem/GO/"))
go_hb <- read_GOs("hb", clusters_hb, paste0(folder, "hb", cl_number[3], "cl_", mem, "mem/GO/"))
go_mb <- read_GOs("mb", clusters_mb, paste0(folder, "mb", cl_number[2], "cl_", mem, "mem/GO/"))

# Extract selected GOs
extract_GO <- function(data, term_pattern, term_description, tissue){
  data_per <- data.frame()
  for (i in 1:length(term_pattern)) {
    data_per1 <- data[data$GOs %in% data$GOs[grep(term_pattern[i], data$GOs)],]
    tryCatch({
      data_per1$Description <- term_description[i]
    }, error=function(e){})
    data_per <- rbind(data_per, data_per1)
  }
  data_per$tissue <- tissue
  data_per
}

subset_go_mb <- extract_GO(go_mb, term_pattern, term_description, "M")
subset_go_hb <- extract_GO(go_hb, term_pattern, term_description, "H")
subset_go_fb <- extract_GO(go_fb, term_pattern, term_description, "F")

# Count occurrences of each term
all_gos <- rbind(rbind(subset_go_fb, subset_go_hb), subset_go_mb)
for (term in term_description) {
  l <- length(unique(all_gos[all_gos$Description == term, 1]))
  all_gos[all_gos$Description == term, 4] <- paste0(term, "\nn=", l)
}

### Change max Benjamini value
all_gos2 <- all_gos
all_gos2$Benjamini[all_gos2$Benjamini > 10] <- 10


order <- names(sort(a))
all_gos2$Condition <- factor(all_gos2$Condition, levels = order)

## Add pattern info
patterns <- a
all_gos2$pattern <- NA
for (i in 1:nrow(all_gos2)) {
  all_gos2$pattern[i] <- patterns[names(patterns) == all_gos2$Condition[i]]
}
all_gos2$pattern <- factor(all_gos2$pattern, levels = 1:max(a))

# plot: dot plot
p <- ggplot(data = all_gos2, aes(x = Condition, y = GOs, color = Benjamini, )) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red", 
                       name = expression(paste("-log"[10],"(FDR)", sep = "")),
                       breaks=c(2.5,5.0,7.5,10),labels=c("2.5","5","7.5", "\u2265 10")) +
  theme_bw() + 
  ylab("") + 
  xlab("") +
  theme(axis.text.y=element_blank(), axis.ticks=element_blank(), panel.grid.major=element_blank(),
        panel.border = element_rect(size = .5, colour = "grey"),
        strip.text.x = element_text(size = 12, face = "bold"),
        strip.text.y = element_text(size = 12),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_x_discrete(drop=TRUE) + scale_fill_discrete(drop=FALSE)

pdf(paste0(folder, "go_tissues.pdf"), width = 8, height = 6) # text might be too big, better plot in R graphs
p + 
  facet_nested(Description ~ pattern + tissue, drop = FALSE, scales = "free", switch = "y") +
  theme(panel.spacing = unit(0, "line"))
dev.off()

# Effects of applying threshold to max padj value
original <- all_gos[, c(1, 3)]
original$distribution <- "original"
new <- all_gos2[, c(1, 3)]
new$distribution <- "thresholded"
df_ggplot <- rbind(original, new)

library(ggrepel)
pdf(paste0(folder, "thresholded_vs_original.pdf"), width = 8, height = 6)
ggplot(df_ggplot, aes(x=Benjamini, fill=distribution)) +
  geom_density(alpha=0.5) +
  labs(y="", x="-log10(FDR)") +
  geom_text_repel(data = subset(df_ggplot, Benjamini>10),  
                  label=subset(df_ggplot, Benjamini>10)[,1],
                  y=0,
                  nudge_y = 0.1,
                  direction = "x",
                  angle        = 90,
                  vjust        = 0,
                  segment.size = 0.2) +
  annotate("text", x = 17.5, y = 0.45, label = paste0("n = ", nrow(df_ggplot)), size=7) +
  theme(panel.background = element_blank(), text = element_text(size=20),
        axis.text = element_text(colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size=.75)) +
  xlim(0, 20)
dev.off()
