# Make two tables:
# first one connecting each significant gene with its associated ATAC peaks
# second one  containing all peaks and stats from those peaks

# lrt version 

tissue <- "FB"
rna <- read.table(paste0("~/nfur/RNAseq/deseq results/lrt/", tissue, "_LRT.txt"), header = T)

atac <- read.table(paste0("~/nfur/ATACseq/edger/diff_peaks_", tissue, ".txt"), header = T)

atac_down_downstream <- read.table(paste0("~/nfur/ATACseq/edger/diff_peaks_down_", tissue, ".downstream.txt"))
atac_down_downstream$coor <- paste0(atac_down_downstream$V1,":", atac_down_downstream$V2, "-", atac_down_downstream$V3)
atac_down_downstream$origin <- "atac_down_downstream"

atac_down_upstream <- read.table(paste0("~/nfur/ATACseq/edger/diff_peaks_down_", tissue, ".upstream.txt"))
atac_down_upstream$coor <- paste0(atac_down_upstream$V1,":", atac_down_upstream$V2, "-", atac_down_upstream$V3)
atac_down_upstream$origin <- "atac_down_upstream"

atac_up_downstream <- read.table(paste0("~/nfur/ATACseq/edger/diff_peaks_up_", tissue, ".downstream.txt"))
atac_up_downstream$coor <- paste0(atac_up_downstream$V1,":", atac_up_downstream$V2, "-", atac_up_downstream$V3)
atac_up_downstream$origin <- "atac_up_downstream"

atac_up_upstream <- read.table(paste0("~/nfur/ATACseq/edger/diff_peaks_up_", tissue, ".upstream.txt"))
atac_up_upstream$coor <- paste0(atac_up_upstream$V1,":", atac_up_upstream$V2, "-", atac_up_upstream$V3)
atac_up_upstream$origin <- "atac_up_upstream"

atac_coor <- rbind(atac_down_downstream, atac_down_upstream, atac_up_downstream, atac_up_upstream)

int_down <- read.table(paste0("~/nfur/intersection_rnaseq_atac/simple/", tissue, "/de.txt"))
int_up <- read.table(paste0("~/nfur/intersection_rnaseq_atac/simple/", tissue, "/ue.txt"))
int <- rbind(int_down, int_up)
int <- rna$gene_id

library(dplyr)
gn_all <- read.table("/home/areyes/nfur/RNAseq/Nfur_gene_names.txt")[,c(1,3)]
gn_unique <- gn_all %>% distinct()
gid <- int
gn <- gn_unique[gn_unique[,1] %in% gid,]
colnames(gn) <- c("gene_id", "V2")

gene_coor_map <- data.frame(gene=gid, coor=NA)
for (i in 1:dim(gene_coor_map)[1]) {
  gene_coor_map[i,2] <- paste(atac_coor[atac_coor$V8 == gene_coor_map[i,1], 10], collapse = "|")
}

table <- data.frame(nfur_name = rna$gene_id,
                    danre_name = merge(gn, rna, by="gene_id")[,2],
                    rna_LFC = rna$log2FoldChange,
                    rna_padj = rna$padj,
                    atac_peak = gene_coor_map$coor)

filt_table <- table[table$atac_peak != "",]

write.table(filt_table, paste0("~/nfur/intersection_rnaseq_atac/simple/lrtdeseq2rna_edgeratac/", tissue, "_diffgene&peaks.tab"), row.names = F, quote = F, sep = "\t")
write.table(atac, paste0("~/nfur/intersection_rnaseq_atac/simple/lrtdeseq2rna_edgeratac/", tissue, "_infopeaks.tab"), row.names = F, quote = F, sep = "\t")