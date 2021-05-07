# SEARCH GENES IN AGING DATABASES

# read databases
cellAge <- read.table("~/nfur/aging_dbs/db/cellAge1.csv", header = T, sep = ";")
agehuman <- read.table("~/nfur/aging_dbs/db/genage_human.csv", sep = ",", header = T)
#agemodels <- read.table("~/nfur/aging_dbs/db/genage_models.csv", sep = ",", header = T)
interactions <- read.delim("~/nfur/aging_dbs/db/interactions.tsv", sep = "\t", header = T)
siginteractions <- interactions[na.omit(interactions$interaction_group_score) > 3,] # Filter score of interaction
drugage <- read.delim("~/nfur/aging_dbs/db/uniquedrugsaging.csv", header = F)[,1]

# Genes lists
geneinteractingagingdrugs <- unique(siginteractions[toupper(drugage) %in% siginteractions$drug_claim_name, 1]) # 983 genes
genesenescence <- cellAge$gene_name # 279 genes
geneshuman <- agehuman$symbol # 307 genes

humanor <- read.table("~/nfur/genome/Nfu_20150522.humanensemblids.tab", header = T)
mygenelist_killi <- read.table("~/nfur/RNAseq/")
mygenelist <- c("ZNF816", "NANOG", "RREB1", "ZBTB4", "NEUROG2", "SREBF2", "ZNF563", "ZFP161", "ZNF331", "RELB", "IRF3", "ZNF232", "ELK1")
mygenelist_subset <- mygenelist[mygenelist %in% humanor$hgnc_symbol]
htest_pval <- function(mygenelist, database){
  background <- humanor$hgnc_symbol
  
  m <- sum(database %in% background)            # effective number background genes
  n <- sum(!(background %in% database))         # number of genes not in background
  k <- length(mygenelist)                       # number of genes in my list of interest
  jointgroup <- sum(mygenelist %in% database)   # number of genes in my list of interest and in the database
  x <- c(0:k)
  
  probabilities <- dhyper(x, m, n, k, log = FALSE)
  pvalue <- sum(probabilities[(jointgroup+1):length(probabilities)])
  return(pvalue)
}

htest_pval(mygenelist = mygenelist_subset, database = geneinteractingagingdrugs)
mygenelist_subset[mygenelist_subset %in% geneinteractingagingdrugs]
htest_pval(mygenelist = mygenelist_subset, database = genesenescence)
mygenelist_subset[mygenelist_subset %in% genesenescence]
htest_pval(mygenelist = mygenelist_subset, database = geneshuman)
mygenelist_subset[mygenelist_subset %in% geneshuman]
