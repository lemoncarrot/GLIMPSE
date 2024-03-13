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

library(googleCloudStorageR)

#Sys.setenv(GCS_AUTH_FILE = "/home/dbr_journalclub/amazing-craft-416704-a9265357ed48.json")
#library(googleCloudStorageR)
#gcs_auth()

#bucket_name <- "writebucket"
#object_name <- "meth_1.png"
#file_path <- "/home/dbr_journalclub/meth_1.png"

#gcs_upload(file_path, bucket = bucket_name, name = object_name)

#convert cpg to genes
library(FDb.InfiniumMethylation.hg19)
filter01 <- read.csv("DMR_0_1.csv")
filter02 <- read.csv("DMR_0_2.csv")
sig.cpg <- filter02$SiteNames
hm450 <- get450k()
probes <- hm450[sig.cpg]
genes <- getNearestTSS(probes)
df <- as.data.frame(genes)
write.csv(df, "DMR_genes_1.csv")

#convert df column to pasteable format (for shinygo)
values <- df$nearestGeneSymbol

values_collapsed <- paste(values, collapse = ",")
writeLines(values_collapsed, "paste2.txt")
