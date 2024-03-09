#converting IDs code
library(dplyr)
library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genes <- getBM(filters="hgnc_symbol", 
               attributes=c("ensembl_gene_id","hgnc_symbol"), 
               values=single_gene_list, 
               mart=ensembl)
mapping_vector <- setNames(genes$ensembl_gene_id, genes$hgnc_symbol)
ensembl_ids <- mapping_vector[single_gene_list]
#ensembl_ids <- as.data.frame(ensembl_ids)
#ensembl_ids <- na.omit(ensembl_ids)
common_ids <- intersect(ensembl_ids, bulk_gene_list)
#write.csv(common_ids, "common_genes.csv", row.names=FALSE)
