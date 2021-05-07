library(biomaRt)
library(dplyr)

writeFasta<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"ensembl_gene_id"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"cdna"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

# set up the database and dataset
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "nfurzeri_gene_ensembl")

# Get all genes from biomart killi assembly 
b <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"), mart=ensembl)
b_seqs <- getSequence(id = b[,1], mart = ensembl, seqType = "cdna", type = "ensembl_gene_id")
writeFasta(b_seqs, "~/nfur/genome/nfur_ensemblid_genes.fasta")

# Get biomart killi genes which are orthologs to human genes
filters <- listFilters(ensembl)
attributes <- listAttributes(ensembl)
attributes[grep("hgnc", attributes$name),] # use this line to explore the list of attributes
attributes_subset <- c("ensembl_gene_id", "hsapiens_homolog_ensembl_gene", "hsapiens_homolog_orthology_type", "hsapiens_homolog_orthology_confidence")
a <- getBM(attributes = attributes_subset, mart=ensembl)
b <- a[a$hsapiens_homolog_ensembl_gene != "",] # 18546 homologs

# Check if retrieving killi homologs from human dataset give same result
ensembl_human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "http://www.ensembl.org")
attributes_subset <- c("ensembl_gene_id", "hgnc_symbol")
a_h <- getBM(attributes = attributes_subset, mart=ensembl_human)
attributes_subset <- c("ensembl_gene_id", "nfurzeri_homolog_ensembl_gene")
c_h <- getBM(attributes = attributes_subset, "ensembl_gene_id", a_h[,1], mart = ensembl_human)
b_h <- merge(a_h, c_h, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id")
b_h <- b_h[b_h$nfurzeri_homolog_ensembl_gene != "" & b_h$hgnc_symbol!= "",] # 18546 homologs not taking into account the hgnc symbol
c <- read.table("~/nfur/genome/Nfu_20150522.ensemblids.tab")
d_h <- merge(b_h, c, by.x="nfurzeri_homolog_ensembl_gene", by.y="V1")
colnames(d_h)[4] <- "own_assembly"
write.table(d_h, "~/nfur/genome/Nfu_20150522.humanensemblids.tab", quote = F, row.names = F, sep = "\t")


