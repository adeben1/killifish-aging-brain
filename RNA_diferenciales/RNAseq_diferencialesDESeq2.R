library(DESeq2)
library(pheatmap)
library(ggplot2)
library(yarrr)
library(PCAtools)
library(RColorBrewer)

folder <- "/home/areyes/nfur/RNAseq/deseq results/wald_test/" # output to be stored
file_counts <- "/home/areyes/nfur/RNAseq/counts/rna_ts_counts_fb.tab"
tissue <- toupper(unlist(strsplit(unlist(strsplit(file_counts, "_"))[4], ".t"))[1])

cts <- as.matrix(read.table(file_counts, header=TRUE, sep="\t", as.is=TRUE, check.names = F))
exclude <- c("__no_feature", "__ambiguous", "__too_low_aQual", "__not_aligned", "__alignment_not_unique")
cts <- cts[!rownames(cts) %in% exclude, ]

coldata <- data.frame(time = rep(c("week8", "week12", "week16", "week20", "week24"), 1, each= 2),
                      rep = rep(c("1", "2"), 5),
                      row.names = colnames(cts),
                      stringsAsFactors = TRUE)

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~rep + time)

dds <- DESeq(dds)
rld <- rlog(dds)
mat <- assay(rld)

# PCA
pdf(paste(folder, "PCA_", tissue, ".pdf", sep=""))
metadata <- data.frame(Age=as.character(c(8, 8, 12, 12, 16, 16, 20, 20, 24, 24)))
rownames(metadata) <- colnames(mat)
p <- pca(mat, metadata = metadata, removeVar = 0.1)
biplot(p, colby = 'Age', legendPosition = "right", lab = NULL, pointSize = 6, legendLabSize = 25, legendIconSize = 8, labSize = 6, axisLabSize = 25) +
  scale_color_brewer(palette = "Set2")
dev.off()

# Plot distribution of expression level
pdf(paste(folder, "Expressionlevels_", tissue, ".pdf", sep=""))
par(mfcol = c(2, 1), mar = numeric(4), oma = c(3, 3, .5, .5), 
    mgp = c(2, .6, 0))
data <- data.frame(expression = as.vector(log2(counts(dds) + 1)),
                  sample = rep(colnames(mat), each=nrow(mat)))
pirateplot(formula = expression ~ sample, data = data, point.o = 0, xaxt = 'n',
           ylab = "", xlab = "", yaxt='n', quant = c(.25, .75),
           quant.col = "black", bean.f.col = "white", theme=3,
           bean.b.col = rep(brewer.pal(n = 5, name = "Dark2"), each=2))
text(9.7, 19.7, "Before normalization", cex = 0.9)
text(1.7, 19.7, paste0("n=", nrow(dds)))
axis(side = 2, las = 2, mgp = c(3, 0.75, 0), yaxp  = c(0, 18, 1))
data <- data.frame(expression = as.vector(log2(counts(dds, normalized=T) + 1)),
                   sample = rep(colnames(mat), each=nrow(mat)))
pirateplot(formula = expression ~ sample, data = data, point.o = 0, xaxt = 'n',
           ylab = "", xlab = "", yaxt='n', quant = c(.25, .75),
           quant.col = "black", bean.f.col = "white", theme=3,
           bean.b.col = rep(brewer.pal(n = 5, name = "Dark2"), each=2))
text(9.7, 19.7, "After normalization", cex = 0.9)
text(.7, 19.7, paste0("n=", nrow(dds)))
axis(side = 2, las = 2, mgp = c(3, 0.75, 0), yaxp  = c(0, 18, 1))
text(x = 1:ncol(mat), y = par("usr")[3] - 1, labels = colnames(mat), xpd = NA,
     srt = 35, adj = 0.8, cex = .9)
mtext(expression(paste("log"[2],"(counts+1)", sep='')), side = 2, outer = TRUE, line = 1.7)
dev.off()

# Extract DEGs
res1 <- results(dds, contrast = c("time", "week12", "week8"))
res2 <- results(dds, contrast = c("time", "week16", "week8"))
res3 <- results(dds, contrast = c("time", "week20", "week8"))
res4 <- results(dds, contrast = c("time", "week24", "week8"))

padj = 0.05
lfc = 0.75

up1 <- res1[res1$log2FoldChange > lfc & res1$padj < padj & !is.na(res1$padj),]
up2 <- res2[res2$log2FoldChange > lfc & res2$padj < padj & !is.na(res2$padj),]
up3 <- res3[res3$log2FoldChange > lfc & res3$padj < padj & !is.na(res3$padj),]
up4 <- res4[res4$log2FoldChange > lfc & res4$padj < padj & !is.na(res4$padj),]

down1 <- res1[res1$log2FoldChange < -lfc & res1$padj < padj & !is.na(res1$padj),]
down2 <- res2[res2$log2FoldChange < -lfc & res2$padj < padj & !is.na(res2$padj),]
down3 <- res3[res3$log2FoldChange < -lfc & res3$padj < padj & !is.na(res3$padj),]
down4 <- res4[res4$log2FoldChange < -lfc & res4$padj < padj & !is.na(res4$padj),]

write.table(rownames(up1), paste(folder, tissue, "_12v8sem_up.txt", sep = ""), quote = F, row.names = F, col.names = F)
write.table(rownames(up2), paste(folder, tissue, "_16v8sem_up.txt", sep = ""), quote = F, row.names = F, col.names = F)
write.table(rownames(up3), paste(folder, tissue, "_20v8sem_up.txt", sep = ""), quote = F, row.names = F, col.names = F)
write.table(rownames(up4), paste(folder, tissue, "_24v8sem_up.txt", sep = ""), quote = F, row.names = F, col.names = F)
write.table(rownames(down1), paste(folder, tissue, "_12v8sem_down.txt", sep = ""), quote = F, row.names = F, col.names = F)
write.table(rownames(down2), paste(folder, tissue, "_16v8sem_down.txt", sep = ""), quote = F, row.names = F, col.names = F)
write.table(rownames(down3), paste(folder, tissue, "_20v8sem_down.txt", sep = ""), quote = F, row.names = F, col.names = F)
write.table(rownames(down4), paste(folder, tissue, "_24v8sem_down.txt", sep = ""), quote = F, row.names = F, col.names = F)

write.table(down4, paste0(folder, tissue, "_24v8sem_down.tab"), quote=F)
write.table(up4, paste0(folder, tissue, "_24v8sem_up.tab"), quote=F)

select_up <- match(rownames(up4), rownames(dds))
select_down <- match(rownames(down4), rownames(dds))
pdf(paste(folder, "Heatmap_", tissue, ".pdf", sep=""))
pheatmap(mat[c(select_up, select_down),], cluster_rows=T, show_rownames=FALSE,
         cluster_cols=FALSE, border_color = NA, scale = "row", clustering_method = "mcquitty", 
         fontsize = 30, labels_col = c("w8_1", "w8_2", "w12_1", "w12_2", "w16_1", "w16_2", "w20_1", "w20_2", "w24_1", "w24_2"))
dev.off()


## Using LRT 
folder <- "/home/areyes/nfur/RNAseq/deseq results/lrt/" # output to be stored
file_counts <- "/home/areyes/nfur/RNAseq/counts/rna_ts_counts_fb.tab"
tissue <- toupper(unlist(strsplit(unlist(strsplit(file_counts, "_"))[4], ".t"))[1])

cts <- as.matrix(read.table(file_counts, header=TRUE, sep="\t", as.is=TRUE, check.names = F))
exclude <- c("__no_feature", "__ambiguous", "__too_low_aQual", "__not_aligned", "__alignment_not_unique")
cts <- cts[!rownames(cts) %in% exclude, ]

coldata <- data.frame(time = rep(c("week8", "week12", "week16", "week20", "week24"), 1, each= 2),
                      rep = rep(c("1", "2"), 5),
                      row.names = colnames(cts),
                      stringsAsFactors = TRUE)

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ rep + time)

# Desing from https://support.bioconductor.org/p/113630/
dds <- DESeq(dds, test = "LRT", reduced = ~rep)
rld <- rlog(dds)
mat <- assay(rld)

res <- results(dds)

padj = 0.001

diff <- res[res$padj < padj & !is.na(res$padj),]

write.table(data.frame("gene_id" = rownames(diff), as.data.frame(diff)), 
            file = paste(folder, tissue, "_LRT.txt", sep = ""), 
            quote = F, sep = "\t", row.names = F)

select_up <- match(rownames(up), rownames(dds))
select_down <- match(rownames(down), rownames(dds))
select <- match(rownames(sig_modelwithrep), rownames(dds))
pdf(paste(folder, "Heatmap_", tissue, ".pdf", sep=""), )
pheatmap(mat[select,], cluster_rows=T, show_rownames=FALSE,
         cluster_cols=FALSE, border_color = NA, scale = "row", clustering_method = "ward.D", 
         fontsize = 30, labels_col = c("w8_1", "w8_2", "w12_1", "w12_2", "w16_1", "w16_2", "w20_1", "w20_2", "w24_1", "w24_2"))
dev.off()

## put (absolute) correlations on the upper panels,
## with size proportional to the correlations.
colnames(mat) <- c("w8_1", "w8_2", "w12_1", "w12_2", "w16_1", "w16_2", "w20_1", "w20_2", "w24_1", "w24_2")
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y, method="spearman"))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
pairs(mat, lower.panel = panel.smooth, upper.panel = panel.cor,
      gap=0, row1attop=FALSE, col= "#00000033", pch=20, lwd=1)

