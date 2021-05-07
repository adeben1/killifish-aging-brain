# Specific genes computed with DESeq2

library(DESeq2)
library(ggplot2)
source("~/nfur/RNAseq/scripts/aux_functions.R")


file_counts <- c("/home/areyes/nfur/RNAseq/counts/rna_ts_counts_fb.tab",
                 "/home/areyes/nfur/RNAseq/counts/rna_ts_counts_hb.tab",
                 "/home/areyes/nfur/RNAseq/counts/rna_ts_counts_mb.tab")

cts_fb <- as.matrix(read.table(file_counts[1], header=TRUE, sep="\t", as.is=TRUE, check.names = F))
cts_mb <- as.matrix(read.table(file_counts[2], header=TRUE, sep="\t", as.is=TRUE, check.names = F))
cts_hb <- as.matrix(read.table(file_counts[3], header=TRUE, sep="\t", as.is=TRUE, check.names = F))

# check datasets are sorted 
fb_mb_order <- all(rownames(cts_fb) == rownames(cts_mb))
fb_hb_order <- all(rownames(cts_fb) == rownames(cts_hb))
if (!fb_mb_order | !fb_hb_order) {
  stop("Can not merge datasets, sort them.")
}

# merge datasets (FB HB MB)
cts <- data.frame()
cts <- cbind(cbind(cts_fb, cts_mb), cts_hb)

# remove rows not containing gene ids
exclude <- rownames(cts)[grep("Nfu_g", rownames(cts), invert = T)]
cts <- cts[!rownames(cts) %in% exclude, ]

# metadata
coldata <- data.frame(age = rep(c(8, 12, 16, 20, 24), 3, each= 2), # age as numeric
                      rep = rep(c("1", "2"), 15),
                      tissue = rep(c("fb", "hb", "mb"), 1, each=10),
                      row.names = colnames(cts),
                      stringsAsFactors = TRUE)

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ tissue + age + tissue:age)

keep <- rowSums(counts(dds) >= 10) >= 20
dds <- dds[keep,]

# run DESeq2
# Genes with small p values from this test are those which at one or more 
# time points after age 8 showed a tissue-specific effect
dds <- DESeq(dds, test="LRT", reduced = ~ tissue + age)

# sometimes the GLM don't converge when there is a single count in a row of all 0's
# to see those genes rownames(dds)[!mcols(dds)$fullBetaConv]
# I choose to omit those genes
ddsClean <- dds[which(mcols(dds)$fullBetaConv),]

res <- results(ddsClean)
head(res[order(res$padj),], 10)
genes_timenumeric <- rownames(res[res$padj<0.05,]) 

gene_plot_deseq2("Nfu_g_1_011301", ddsClean)

fiss <- plotCounts(dds, "Nfu_g_1_011301", 
                   intgroup = c("age","tissue"), returnData = TRUE)
fiss$age <- as.numeric(as.character(fiss$age))
ggplot(fiss,
       aes(x = age, y = count, color = tissue, group = tissue)) + 
  geom_point() + stat_summary(fun=mean, geom="line") +
  scale_y_log10() + ggtitle("CDH6") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=25)) +
  labs(y = "Normalized DESeq2 counts", x = "Age (weeks)") +
  scale_x_continuous(breaks = c(8, 12, 16, 20, 24)) +
  scale_color_manual(values=c("green4", "orange1", "blue"))
