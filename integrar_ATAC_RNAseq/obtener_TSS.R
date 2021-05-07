# Previosly:
#           grep -P "\tgene\t" Nfu_20150522.genes_20150922.gff3 | cut -f1,4,5,7 > Nfu_20150522.genes_20150922.TSS.tmp.bed
#           grep -P "\tgene\t" Nfu_20150522.genes_20150922.gff3 | cut -f9 > Nfu_20150522.genes_20150922.TSS.tmp.descriptions
#           cut -d ';' -f1 Nfu_20150522.genes_20150922.TSS.tmp.descriptions | cut -d '=' -f2 > Nfu_20150522.genes_20150922.TSS.names
#           paste Nfu_20150522.genes_20150922.TSS.tmp.bed Nfu_20150522.genes_20150922.TSS.names > final

setwd("nfur/")

gtf <- read.table("Nfu_20150522.genes_20150922.TSS.tmp.bed", sep = "\t")
window <- 50 # 50*2=100 bps

for (line in seq(1:nrow(gtf))) {
  if (gtf[line, 4] == "-") {
    tss <- gtf[line, 3]
  }
  else {
    tss <- gtf[line, 2]
  }
  
  gtf[line, 3] <- tss + window
  if (tss < 50) {
    gtf[line, 2] <- 0
  }
  else {
    gtf[line, 2] <- tss - window
  }
}

write.table(gtf, "Nfu_20150522.genes_20150922.TSS.bed", sep = "\t", quote = F, col.names = F, row.names = F)

# sortBed -i Nfu_20150522.genes_20150922.TSS.bed > Nfu_20150522.genes_20150922.TSS.sorted.bed