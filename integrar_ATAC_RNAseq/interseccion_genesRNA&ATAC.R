##### INTERSECT RNASEQ DESEQ2 AND ATASEQ EDGER ###### 
library(eulerr)

rna_folder <- "/home/areyes/nfur/RNAseq/deseq results/lrt/"
atac_folder <- "/home/areyes/nfur/ATACseq/edger/"

hb <- list(RNAseq = unique(read.table(paste0(rna_folder, "HB_LRT.txt"), header = T)[,1]),
           ATACseq_down = unique(read.table(paste0(atac_folder, "diff_peaks_down_HB.killigenes.txt"))[,1]),
           ATACseq_up = unique(read.table(paste0(atac_folder, "diff_peaks_up_HB.killigenes.txt"))[,1]))
hb$ATACseq <- unique(c(hb$ATACseq_down, hb$ATACseq_up))
plot(euler(hb[c(1,4)]), quantities=T, main = "Hind brain")

mb <- list(RNAseq = unique(read.table(paste0(rna_folder, "MB_LRT.txt"), header = T)[,1]),
           ATACseq_down = unique(read.table(paste0(atac_folder, "diff_peaks_down_MB.killigenes.txt"))[,1]),
           ATACseq_up = unique(read.table(paste0(atac_folder, "diff_peaks_up_MB.killigenes.txt"))[,1]))
mb$ATACseq <- unique(c(mb$ATACseq_down, mb$ATACseq_up))
plot(euler(mb[c(1,4)]), quantities=T, main = "Mid brain")

fb <- list(RNAseq = unique(read.table(paste0(rna_folder, "FB_LRT.txt"), header = T)[,1]),
           ATACseq_down = unique(read.table(paste0(atac_folder, "diff_peaks_down_FB.killigenes.txt"))[,1]),
           ATACseq_up = unique(read.table(paste0(atac_folder, "diff_peaks_up_FB.killigenes.txt"))[,1]))
fb$ATACseq <- unique(c(fb$ATACseq_down, fb$ATACseq_up))
plot(euler(fb[c(1,4)]), quantities=T, main = "Fore brain")

# Save files
source("/home/areyes/nfur/RNAseq/scripts/aux_functions.R")

intersect_genes <- function(tissue, name){
  intersect <- intersect(tissue$RNAseq, tissue$ATACseq)

  dir_name <- paste0("~/nfur/intersection_rnaseq_atac/simple/", name)
  dir.create(dir_name)
  setwd(dir_name)
  
  write.table(intersect, "intersect.txt", sep="\t", col.names = F, row.names = F, quote = F)
  
  nfugene_to_dangene(dir_name)
  
  pantherGO(dir_name)
}

intersect_genes(fb, "FB_lrt_edger")
intersect_genes(hb, "HB_lrt_edger")
intersect_genes(mb, "MB_lrt_edger")


## GO PLOT
filename <- "/home/areyes/nfur/intersection_rnaseq_atac/simple/MB_lrt_edger/GO/intersect_GO.txt"

go_table <- read.table(filename, header = T, sep = "\t")
go_table <- go_table[!go_table$term.label %in% c("cellular process", "UNCLASSIFIED", "biological_process"),]
go_table$fdr <- -log10(go_table$fdr)

# In case some term is too large
new_term <- "regulation... metabolic process"
go_table$term.label <- unfactor(go_table$term.label)
go_table$term.label[go_table$term.label %in% "regulation of nucleobase-containing compound metabolic process"] <- new_term

library(RColorBrewer)
color <- brewer.pal(n = 5, name = "Set2")

library(ggplot2)
ggplot(head(go_table), aes(x=reorder(term.label, -fdr), y=fdr)) +
  geom_col(width = 0.7, show.legend = F, fill=color[5]) +
  coord_flip() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.line.x.bottom = element_line(),
        axis.line.y.left = element_line(),
        axis.title.y = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        text = element_text(size=17), 
        aspect.ratio = 1, 
        axis.text.y = element_text(size=15)) +
  labs(title= "MB ATAC&RNA", y = expression("-log"[10]* "(FDR)"))

ggsave(paste0(filename, ".pdf"), height = 6, width = 6)
