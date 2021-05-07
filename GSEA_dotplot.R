# Plot results from GSEA

library(ggplot2)
library(dplyr)

root.folder <- "~/nfur/RNAseq/GSEA/"
gsea_result_folders <- list.files(root.folder, pattern = glob2rx("deseq2*hallmarks*"))
gsea_result_files <- lapply(gsea_result_folders, function(x) list.files(paste0(root.folder, x), pattern = glob2rx("*report*tsv"), full.names = T))

# Read NAME, NES and FDR
results_raw <- lapply(unlist(gsea_result_files), function(x) read.table(file = x, sep = "\t", header = T)[,c(1, 6, 8)])
results_fil <- lapply(results_raw, function(x) x[x$FDR.q.val<=0.05,])

# Add tissue data
results_fil[[1]] <- cbind(results_fil[[1]], tissue = "Fore brain")
results_fil[[2]] <- cbind(results_fil[[2]], tissue = "Fore brain")
results_fil[[3]] <- cbind(results_fil[[3]], tissue = "Hind brain")
results_fil[[4]] <- cbind(results_fil[[4]], tissue = "Hind brain")
results_fil[[5]] <- cbind(results_fil[[5]], tissue = "Mid brain")
results_fil[[6]] <- cbind(results_fil[[6]], tissue = "Mid brain")

# Merge in one table
pathway_table <- bind_rows(results_fil)

# Sort by NES using factors
pathway_table$tissue <- factor(pathway_table$tissue, levels = c("Fore brain", "Mid brain", "Hind brain"))

pathway_table$NAME <- gsub("HALLMARK_", "", pathway_table$NAME)
pathway_table$NAME <- gsub("_", " ", pathway_table$NAME)

# Dotplot
ggplot(pathway_table, aes(x=tissue, y=reorder(NAME, NES))) +
  geom_point(aes(color=NES, size=FDR.q.val)) +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size("FDR", trans = "reverse", limits = c(0.05, 0)) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1))
