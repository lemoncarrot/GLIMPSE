library(Seurat)
library(dplyr)
library(DEsingle)
library(DESeq2)
library(SingleR)
library(celldex)
library(edgeR)
library(ggplot2)
library(pheatmap)
library(InteractiveComplexHeatmap)

mrna_693_normal <- readRDS("mrna_693_normal.rds")
mrna_693_normal$cell_type <- "normal"
testing <- read.csv("mrna_693_final2.csv", header=TRUE)
testing <- testing[, -1]
testing$cell_type <- "tumor"
intersection <- Reduce(intersect, list(colnames(mrna_693_normal), colnames(testing)))
mrna_693_normal <- subset(mrna_693_normal, select=intersection)
testing <- subset(testing, select=intersection)
combined_df <- rbind(testing, mrna_693_normal)
contains_ENSG <- grepl("ENSG", names(combined_df))
cts <- as.matrix(t(combined_df[, contains_ENSG]))
row_names <- rownames(cts)
col_names <- colnames(cts)
cts <- matrix(as.numeric(cts), nrow=nrow(cts), ncol=ncol(cts), dimnames = list(row_names, col_names))
coldata <- data.frame(cell_type = combined_df$cell_type)
rownames(coldata) <- colnames(cts)
#coldata format: 1 column titled "celltype", with cell types listed (matches to columns in cts)
#cts format: genes as rows, columns as samples (in this case, no IDs for samples)

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design= ~ cell_type)

dds <- DESeq(dds)
#interactivate(dds) interactive heatmap
resultsNames(dds) # lists the coefficients
res <- results(dds, name="cell_type_tumor_vs_normal")
#summary(res)
#sum(res$padj < 0.1, na.rm=TRUE)

#MA plot might not be necessary (but scRNA can't get volcano plot)
#png(filename="bulk_1.png", width=8, height=6, units="in", res=300)
#ggplot(res, aes(x=baseMean, y=log2FoldChange, color=(padj<0.05 & abs(log2FoldChange)>2))) + 
#  geom_point(alpha=0.5) + 
#  scale_color_manual(values=c("black", "red")) +
#  labs(title="MA Plot", x="Average expression (baseMean)", y="Log2 fold change") +
#  theme_minimal()
#dev.off()
png(filename="bulk_2.png", width=8, height=6, units="in", res=300)
ggplot(res, aes(x=log2FoldChange, y=-log10(padj), color=(padj<0.05 & abs(log2FoldChange)>2))) + 
  geom_point(alpha=0.5) +
  scale_color_manual(values=c("black", "red")) +
  labs(title=" ", x="Log2 fold change", y="-Log10 adjusted p-value") + 
  theme_minimal() +
  theme(legend.title=element_blank())
dev.off()

sig_genes <- res[!is.na(res$padj) & res$padj < 0.05,]
sig_genes <- sig_genes[abs(sig_genes$log2FoldChange) > 2,]
sig_gene_names <- rownames(sig_genes)
sig_gene_df <- data.frame(Gene = sig_gene_names)
#write.csv(sig_gene_df, "bulk_both_genes.csv", row.names = FALSE, quote = FALSE)
#3188 genes, ENSEMBL IDs (ENSG)
