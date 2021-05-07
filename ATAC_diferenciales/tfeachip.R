library(TFEA.ChIP)

extract_factors_gimme_maelstrom <- function(folder){
  # Input: folder path containing gimme maelstrom results
  
  # read diff motifs
  diff_motif_table <- read.table(paste0(folder, "final.out.txt"), sep = "\t", header = T)
  diff_motif <- diff_motif_table[abs(diff_motif_table[,2]) > 3 |
                                 abs(diff_motif_table[,3]) > 3 |
                                 abs(diff_motif_table[,4]) > 3, 1]
  
  # read gimme database
  gimme_database <- read.table(paste0(folder, "gimme.vertebrate.v5.0.motif2factors.txt"), sep = "\t", header = T)

  # extract factors 
  subset <- data.frame()
  for (motif in diff_motif) {
    subset <- rbind(subset, gimme_database[gimme_database$Motif == motif,])
  }
  
  # keep only curated ones
  subset$Factor <- toupper(subset$Factor)
  curated_factors <- unique(subset[subset$Curated == "Y",2])
  return(curated_factors)
}

diff_expr_raw <- read.table("~/nfur/RNAseq/deseq/HB_LRT.tab", header = T)
human_genes <- read.table("~/nfur/genome/Nfu_20150522.humanensemblids.tab", header = T)[,c(3,4)]

diff_expr_homologs <- merge(diff_expr_raw, human_genes, by.x="genes", by.y="own_assembly")
diff_expr_homologs$genes <- diff_expr_homologs$hgnc_symbol
diff_expr_homologs <- diff_expr_homologs[,-c(3,5)]
names(diff_expr_homologs) <- c("Genes", "log2FoldChange", "pval.adj")

diff_expr_table <- preprocessInputData(diff_expr_homologs)
diff_expr_table <- diff_expr_table[!is.na(diff_expr_table$pval.adj),]

#extract vector with names of upregulated genes
Genes.Upreg <- Select_genes( diff_expr_table, max_pval = 0.001 )
#extract vector with names of non-responsive genes
Genes.Control <- Select_genes( diff_expr_table,
                               min_pval = 0.5, max_pval = 1)

factors <- extract_factors_gimme_maelstrom("~/nfur/ATACseq/gimme_maelstrom/hb_diffpeaks_merged/")
chip_index <- get_chip_index(encodeFilter = T)
CM_list_UP <- contingency_matrix( Genes.Upreg, Genes.Control, chip_index)
pval_mat_UP <- getCMstats( CM_list_UP, chip_index)
pval_mat_UP[pval_mat_UP$adj.p.value<0.05,]

factors[which(factors %in% pval_mat_UP[pval_mat_UP$adj.p.value<0.05,4])]

#TF_ranking <- rankTFs( pval_mat_UP, rankMethod = "gsea", makePlot = TRUE )
#TF_ranking[[ "TFranking_plot" ]]

#plot_CM(pval_mat_UP)
#curated_factors[curated_factors$Factor == "REST",]

