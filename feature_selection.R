library(Seurat)
library(dplyr)
library(DEsingle)
library(DESeq2)
library(SingleR)
library(celldex)
library(edgeR)

x <- readRDS("GSE84465.rds") #scrna dataset
gene_data <- as.data.frame(x[60:nrow(x), ])
metadata <- x[1:59, ]
metadata <- t(metadata)
metadata <- as.data.frame(metadata)
seurat_obj <- CreateSeuratObject(counts = gene_data, project = "GBM-scRNA", meta.data = metadata)

tumor_obj <- subset(seurat_obj, subset = characteristics_ch1.3 == "tissue: Tumor")
peripheral_obj <- subset(seurat_obj, subset = characteristics_ch1.3 == "tissue: Periphery")

# Process Tumor Data
tumor_obj <- subset(tumor_obj, min.cells = 3, min.features=100)
tumor_obj[["percent.mt"]] <- PercentageFeatureSet(tumor_obj, pattern = "^MT-")
tumor_obj <- subset(tumor_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)
#ncol(tumor_obj)
#pdf("scRNA_1.pdf", width=8, height=6)
#VlnPlot(tumor_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#dev.off()
tumor_obj <- NormalizeData(tumor_obj, normalization.method = "LogNormalize", scale.factor = 10000)
tumor_obj <- FindVariableFeatures(tumor_obj, selection.method = "vst", nfeatures = 1500)
#pdf("scRNA_2.pdf", width = 8, height = 6)
#VariableFeaturePlot(tumor_obj)
#dev.off()
tumor_obj <- ScaleData(tumor_obj)
tumor_obj <- RunPCA(tumor_obj, features = VariableFeatures(object = tumor_obj))
#ElbowPlot(tumor_obj) not needed
tumor_obj <- RunTSNE(tumor_obj, dims = 1:20)
tumor_obj <- FindNeighbors(tumor_obj, dims = 1:20)
tumor_obj <- FindClusters(tumor_obj, resolution = 0.5) #can change resolution to adjust number of clusters
data_for_singleR <- GetAssayData(tumor_obj, layer = "RNA", slot = "data")
reference <- HumanPrimaryCellAtlasData()
annotations <- SingleR(test = data_for_singleR, ref = reference, labels = reference$label.main)
head(annotations)
tumor_obj$SingleR.celltype <- annotations$labels
threshold <- 15 #threshold to remove cell annotations with <15 cells 
cell_type_counts <- table(tumor_obj$SingleR.celltype)
rare_cell_types <- names(cell_type_counts[cell_type_counts < threshold])
tumor_obj$SingleR.celltype_adjusted <- tumor_obj$SingleR.celltype
tumor_obj$SingleR.celltype_adjusted[tumor_obj$SingleR.celltype %in% rare_cell_types] <- 'Other'
#pdf("scRNA_3.pdf", width=8, height=6)
#DimPlot(tumor_obj, reduction = "tsne", group.by = "SingleR.celltype_adjusted", label = TRUE)
#dev.off()
tumor_markers <- FindAllMarkers(tumor_obj, min.pct = 0.25, logfc.threshold = 0.25)
significant_tumor_markers <- subset(tumor_markers, p_val_adj < 0.05)
top10_markers <- significant_tumor_markers %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)
marker_genes <- top10_markers$gene
pdf("scRNA_4.pdf", width = 16, height = 10)
DoHeatmap(tumor_obj, features = marker_genes)
dev.off()
unique_genes <- unique(significant_tumor_markers$gene)
genes_df <- data.frame(Gene = unique_genes)
write.csv(genes_df, "scRNA_t_genes.csv", row.names = FALSE)

# Process Peripheral Data
peripheral_obj <- subset(peripheral_obj, min.cells = 3, min.features=100)
peripheral_obj[["percent.mt"]] <- PercentageFeatureSet(peripheral_obj, pattern = "^MT-")
peripheral_obj <- subset(peripheral_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)
#ncol(peripheral_obj)
pdf("scRNA_5.pdf", width = 8, height = 6)
VlnPlot(peripheral_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
peripheral_obj <- NormalizeData(peripheral_obj, normalization.method = "LogNormalize", scale.factor = 10000)
peripheral_obj <- FindVariableFeatures(peripheral_obj, selection.method = "vst", nfeatures = 1500)
pdf("scRNA_6.pdf", width = 8, height = 6)
VariableFeaturePlot(peripheral_obj)
dev.off()
peripheral_obj <- ScaleData(peripheral_obj)
peripheral_obj <- RunPCA(peripheral_obj, features = VariableFeatures(object = peripheral_obj))
#ElbowPlot(tumor_obj) not needed
peripheral_obj <- RunTSNE(peripheral_obj, dims = 1:20)
peripheral_obj <- FindNeighbors(peripheral_obj, dims = 1:20)
peripheral_obj <- FindClusters(peripheral_obj, resolution = 0.5)
data_for_singleR <- GetAssayData(peripheral_obj, layer = "RNA", slot = "data")
reference <- HumanPrimaryCellAtlasData()
annotations <- SingleR(test = data_for_singleR, ref = reference, labels = reference$label.main)
head(annotations)
peripheral_obj$SingleR.celltype <- annotations$labels
threshold <- 15 #threshold to remove cell annotations with <15 cells 
cell_type_counts <- table(peripheral_obj$SingleR.celltype)
rare_cell_types <- names(cell_type_counts[cell_type_counts < threshold])
peripheral_obj$SingleR.celltype_adjusted <- peripheral_obj$SingleR.celltype
peripheral_obj$SingleR.celltype_adjusted[peripheral_obj$SingleR.celltype %in% rare_cell_types] <- 'Other'
pdf("scRNA_7.pdf", width = 8, height = 6)
DimPlot(peripheral_obj, reduction = "tsne", label = TRUE) #group.by= #change to one with annotations
dev.off()
peripheral_markers <- FindAllMarkers(peripheral_obj, min.pct = 0.25, logfc.threshold = 0.25)
significant_peripheral_markers <- subset(peripheral_markers, p_val_adj < 0.05)
top10_markers2 <- significant_peripheral_markers %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)
marker_genes2 <- top10_markers2$gene
pdf("scRNA_8.pdf", width = 16, height = 10)
DoHeatmap(peripheral_obj, features = marker_genes)
dev.off()
unique_genes <- unique(significant_peripheral_markers$gene)
genes_df <- data.frame(Gene = unique_genes)
write.csv(genes_df, "scRNA_p_genes.csv", row.names = FALSE)

combined_obj <- merge(peripheral_obj, y = tumor_obj, add.cell.ids = c("Peripheral", "Tumor"), project = "CombinedProject")
rna_assay <- combined_obj[["RNA"]]
peripheral_counts <- rna_assay@layers$counts.1
tumor_counts <- rna_assay@layers$counts.2
rownames(peripheral_counts) <- rownames(rna_assay)
rownames(tumor_counts) <- rownames(rna_assay)
combined_counts <- cbind(peripheral_counts, tumor_counts)
counts <- as.matrix(combined_counts)
group <- factor(combined_obj@meta.data$characteristics_ch1.3)

valid_columns <- colSums(counts) > 0
counts_filtered <- counts[, valid_columns]
group_filtered <- group[valid_columns]
dge <- DGEList(counts = counts_filtered)
dge$samples$group <- factor(group_filtered)
keep <- filterByExpr(dge)
dge <- dge[keep,]
dge <- calcNormFactors(dge)
design <- model.matrix(~ 0 + dge$samples$group)
dge$design <- design
dge <- estimateGLMCommonDisp(dge, dge$design)
dge <- estimateGLMTrendedDisp(dge, dge$design)
dge <- estimateGLMTagwiseDisp(dge, dge$design)
fit <- glmFit(dge, dge$design)
lrt <- glmLRT(fit)
top_genes <- topTags(lrt, n=20000)
print(top_genes)

#visualization
#not tested yet
#will be scRNA_9.png
sig_threshold <- 0.05
logFC_threshold <- 1 
sig_genes <- decideTestsDGE(lrt, adjust.method="BH", p.value=sig_threshold)
sig_gene_ids <- which(sig_genes != 0) 
sig_gene_results <- topTags(lrt, n=Inf)$table[sig_gene_ids, ]
library(ggplot2)
ggplot(as.data.frame(lrt$table), aes(x=logFC, y=-log10(PValue))) +
  geom_point(aes(color=adj.P.Val < sig_threshold & abs(logFC) > logFC_threshold), alpha=0.4) +
  scale_color_manual(values=c("grey", "red")) +
  labs(title="Volcano plot of differential expression",
       x="Log2 fold change",
       y="-Log10 p-value") +
  theme_minimal()

#heatmap
#not tested yet
sig_gene_names <- rownames(sig_gene_results)
heatmap_data <- counts_filtered[sig_gene_names, ]
heatmap_data_log <- log2(heatmap_data + 1)
library(pheatmap)
png(filename = "scRNA_10.png", width = 1200, height = 1200, units = "px", res = 300)
pheatmap(heatmap_data_log,
         show_rownames = TRUE,
         show_colnames = TRUE,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         scale = "row", # Scale rows (genes) to have zero mean and unit variance
         annotation_col = data.frame(Group = group_filtered), # Add sample grouping as column annotation
         cluster_cols = FALSE) # Do not cluster columns (samples) as per your requirement
dev.off()






# countData: a matrix of RNA-seq counts, genes as rows and samples as columns
# condition: a factor or character vector indicating the condition of each sample (e.g., "tumor", "normal")
#genes with log2FC>2 and adjusted p-value<0.05
results <- DEsingle(combined_counts, group = factor(cell_labels), parallel=TRUE)
saveRDS(results, "scRNA_results.rds", compress=FALSE)
filtered_results <- subset(results, abs(norm_foldChange) <= 400)
significantGenes <- subset(filtered_results, pvalue.adj.FDR < 0.05 & abs(norm_foldChange) > 2 & abs(norm_foldChange) < 5000)
genes_df <- data.frame(Gene = significantGenes$gene)
write.csv(genes_df, "scRNA_both_genes.csv", row.names = FALSE, quote = FALSE)

gene_number = 20 #change how many significant genes we want to examine
top_de_genes <- head(significantGenes$gene, gene_number)
data_for_heatmap <- FetchData(combined_obj, vars = top_de_genes)
library(pheatmap)
pdf("scRNA_9.pdf", width = 10, height = 8)
pheatmap(t(data_for_heatmap), show_rownames = TRUE, show_colnames = FALSE, clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", clustering_method = "complete", scale = "row", color = colorRampPalette(c("blue", "white", "red"))(255))
dev.off()



#__---------------------------------------------------------------------------------------

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

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design= ~ cell_type)

dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds, name="cell_type_tumor_vs_normal")
sig_genes <- res[!is.na(res$padj) & res$padj < 0.05,]
sig_genes <- sig_genes[abs(sig_genes$log2FoldChange) > 1,]
plot(-log10(res$pvalue), res$log2FoldChange, pch=20, main="Volcano Plot", xlab="-log10 p-value", ylab="Log2 Fold Change")
with(sig_genes, points(-log10(pvalue), log2FoldChange, col="red", pch=20))
save(sig_genes, file = "sig_genes_bulk_obj.RData")
sig_gene_names <- rownames(sig_genes)
sig_gene_df <- data.frame(Gene = sig_gene_names)
write.csv(sig_gene_df, "bulk_both_genes.csv", row.names = FALSE, quote = FALSE)
sig_gene_exp <- cts[sig_gene_names, ]
ordered_genes <- res[order(res$pvalue, decreasing = FALSE), ]
top_genes <- head(ordered_genes, 20000)
top_gene_names <- rownames(top_genes)
top_gene_exp <- cts[top_gene_names, ]
if (!is.null(top_gene_exp) && nrow(top_gene_exp) > 0) {
  norm_exp_data <- log2(top_gene_exp + 1)
  pdf("Top_20000_Genes_Heatmap.pdf", width = 12, height = 10)
  pheatmap(norm_exp_data,
           show_rownames = TRUE,
           show_colnames = FALSE,
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           clustering_method = "complete",
           scale = "row",
           cluster_cols = FALSE,
           color = colorRampPalette(c("blue", "white", "red"))(255))
  dev.off()
} else {
  print("No significant genes for heatmap.\n")
}

if (!is.null(sig_gene_exp) && nrow(sig_gene_exp) > 0) {
  norm_exp_data <- log2(sig_gene_exp + 1)
  pdf("bulk_2.pdf", width = 12, height = 10)
  pheatmap(norm_exp_data,
           show_rownames = TRUE,
           show_colnames = FALSE,
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           clustering_method = "complete",
           scale = "row",
           color = colorRampPalette(c("blue", "white", "red"))(255))
  dev.off()
} else {
  print("No significant genes for heatmap.\n")
}




#plotMA(res, main="MA Plot", ylim=c(-2,2))
#points(sig_genes$log2FoldChange, sig_genes$padj, col="red", pch=20)
#library(pheatmap)
#select_genes <- rownames(sig_genes)
#normalized_counts <- assay(rlogTransformation(dds))
#pheatmap(normalized_counts[select_genes, ], scale="row", clustering_distance_rows="euclidean", clustering_distance_cols="euclidean")

x <- readRDS("DEsingleresults.rds")
x <- x[x$norm_foldChange >= 2 | x$norm_foldChange <= -2, ]
x <- x[x$pvalue.adj.FDR<0.05, ]
x$gene_names <- rownames(x)
single_gene_list <- x$gene_names
load("sig_genes_bulk_obj.RData")
gene_table <- sig_genes@listData
gene_table <- as.data.frame(gene_table)
gene_table <- gene_table[gene_table$log2FoldChange>=2 | gene_table$log2FoldChange <=-2,]
gene_table <- gene_table[gene_table$padj<0.05,]
bulk_gene_list <- gene_table$gene_name
#write.csv(bulk_gene_list, "bulk_gene_list.csv", row.names=FALSE)
#write.csv(single_gene_list, "single_gene_list.csv", row.names=FALSE)
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

library(org.Hs.eg.db)
library(AnnotationDbi)
library(enrichplot)
library(ggplot2)

entrezIDs <- bitr(common_ids, fromType = "ENSEMBLID", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
missingGenes <- setdiff(geneList, entrezIDs$SYMBOL)

goResults <- enrichGO(entrezIDs$ENTREZID, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "ALL")
keggResults <- enrichKEGG(entrezIDs$ENTREZID, organism = "hsa", keyType = "ENTREZID")

barplot(goResults, showCategory=20)  # showing top 20 categories

dotplot(keggResults)

emapplot(goResults)
emapplot(keggResults)


#-------------------------------------------------------------------------------------------------------------------------


#load methylation glioblastoma data
#filter down to GBM and control

#perform DE methylation analysis 

#do GO and KEGG pathway analysis (WebGestalt)
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "AnnotationDbi", "enrichplot", "ggplot2"))
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(enrichplot)
library(ggplot2)

#gene list here

#if gene Ids and not entrezIDs, convert
entrezIDs <- bitr(geneList, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
missingGenes <- setdiff(geneList, entrezIDs$SYMBOL)

goResults <- enrichGO(entrezIDs$ENTREZID, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "ALL")
keggResults <- enrichKEGG(entrezIDs$ENTREZID, organism = "hsa", keyType = "ENTREZID")

barplot(goResults, showCategory=20)  # showing top 20 categories

dotplot(keggResults)

emapplot(goResults)
emapplot(keggResults)

#load dataset with prognosis + methylation + expression profiling
#TCGA-GBM
TCGA-GBM <- readRDS("TCGA-GBM.rds")
#build model here

# Load necessary libraries
library(readr)
library(limma)
library(glmnet)
library(survival)
library(caret)
library(rms)
library(Hmisc)
library(ranger)
library(ggplot2)
library(dplyr)
library(ggfortify)
library(data.table)
library(parallel)
library(doParallel)
library(foreach)
library(survminer) # for ggsurvplot
library(timeROC)


# Load datasets
sig_genes <- read.csv('common_genes.csv')
tcga_rna <- readRDS('TCGA-GBM_rna.rds')
tcga_meth <- readRDS('TCGA-GBM_meth.rds')

# Extract significant gene IDs
gene_ids <- sig_genes$x
base_gene_ids <- sapply(strsplit(gene_ids, "\\."), `[`, 1)

base_gene_ids <- sapply(strsplit(base_gene_ids, "\\."), `[`, 1) # If already done, ensure it matches the RNA-seq data format
colnames(tcga_rna) <- sapply(strsplit(colnames(tcga_rna), "\\."), `[`, 1)
filtered_rna_data <- tcga_rna[, c("sample", "days_to_death", "vital_status", grep(paste0("^", paste(base_gene_ids, collapse="|"), "$"), colnames(tcga_rna), value = TRUE))]

# Filter RNA-seq and methylation data to include relevant columns
filtered_rna_data <- tcga_rna[, c("sample", "days_to_death", "vital_status", base_gene_ids[base_gene_ids %in% colnames(tcga_rna)]), drop = FALSE]
filtered_meth_data <- tcga_meth[, c("sample", "days_to_death", "vital_status", grep("^cg", colnames(tcga_meth), value = TRUE))]


# Identify and retain only overlapping samples
overlap_samples <- intersect(filtered_rna_data$sample, filtered_meth_data$sample)
filtered_rna_data <- filtered_rna_data[filtered_rna_data$sample %in% overlap_samples, ]
filtered_meth_data <- filtered_meth_data[filtered_meth_data$sample %in% overlap_samples, ]

# Ensure no missing values in survival data
complete_cases <- complete.cases(filtered_rna_data$days_to_death, filtered_rna_data$vital_status)
filtered_rna_data <- filtered_rna_data[complete_cases, ]
filtered_meth_data <- filtered_meth_data[complete_cases, ]

combined_data_filtered <- combined_data[, colSums(is.na(combined_data)) < nrow(combined_data)]

threshold_percentage <- 75
# Calculate the percentage of non-NA values for each column
percentage_non_na_per_column <- colSums(!is.na(combined_data_filtered)) / nrow(combined_data_filtered) * 100
# Filter columns based on the threshold
combined_data_filtered_cols <- combined_data_filtered[, percentage_non_na_per_column >= threshold_percentage]

#-----------------
# Identify columns with 100% non-NA values
complete_cols <- colSums(is.na(combined_data_filtered_cols)) == 0

# Subset the data to keep only complete columns
data_complete_cols <- combined_data_filtered_cols[, complete_cols]

# Display the structure of the filtered dataset
str(data_complete_cols)
setDT(data_complete_cols)

# Extract survival-related data into a separate dataframe
survival_data <- data_complete_cols[, .(days_to_death, vital_status)]

# View the first few rows to confirm
head(survival_data)

# Remove metadata columns to prepare for Elastic Net model
feature_data <- data_complete_cols[, !names(data_complete_cols) %in% c("sample", "days_to_death", "vital_status"), with = FALSE]

# View the structure to confirm
str(feature_data)
days_to_death <- data_complete_cols$days_to_death

#------------------
# Combine RNA and Methylation data
rna_features <- setdiff(names(filtered_rna_data), c("sample", "days_to_death", "vital_status"))
meth_features <- setdiff(names(filtered_meth_data), c("sample", "days_to_death", "vital_status"))
X_rna_final <- as.matrix(filtered_rna_data[, rna_features])
X_meth_final <- as.matrix(filtered_meth_data[, meth_features])
X_combined <- cbind(X_rna_final, X_meth_final)

# Order both datasets by 'sample' column to ensure consistent sample alignment
filtered_rna_data <- filtered_rna_data[order(filtered_rna_data$sample), ]
filtered_meth_data <- filtered_meth_data[order(filtered_meth_data$sample), ]

# Exclude the metadata columns from the methylation dataset before merging
# Assuming 'sample' is the first column and you want to keep it from the RNA dataset
meth_features_only <- filtered_meth_data[, -c(1, which(names(filtered_meth_data) %in% c("days_to_death", "vital_status")))]

# Combine the RNA data (with metadata) and methylation features
combined_data <- cbind(filtered_rna_data, meth_features_only)

#---------------
# Convert 'vital_status' into a binary format where, for example, "Dead" is 1 (event occurred) and others are 0 (censored)
vital_status_binary <- as.integer(survival_data$vital_status == "Dead")

# Create the survival object
surv_obj <- Surv(time = survival_data$days_to_death, event = vital_status_binary)

# Now, create a matrix y with 'time' and 'status' columns from surv_obj
# Note: This step assumes glmnet version compatibility, where a matrix format is acceptable
# Extract time and status from surv_obj (assuming surv_obj can be treated like a dataframe for extraction, which it cannot directly. This is for conceptual understanding)
# Assuming filtered_rna_data$days_to_death and vital_status_binary are correctly prepared
time <- survival_data$days_to_death  # Time to event or censoring
status <- vital_status_binary  # Event occurred (1) or censored (0)
y <- matrix(c(time, status), ncol = 2, byrow = FALSE)
colnames(y) <- c("time", "status")

x <- as.matrix(feature_data)


# Assuming 'X_combined' as your feature matrix and 'days_to_death' as your response variable

# Load necessary libraries
library(caret)
library(glmnet)

# Set up cross-validation
control <- trainControl(method = "cv", number = 10, search = "grid")

# Define the grid for hyperparameter tuning
grid <- expand.grid(
  alpha = seq(0, 1, length.out = 5), # Adjust the sequence as needed
  lambda = 10^seq(-3, 3, length.out = 100) # Adjust the sequence as needed
)

# Train the model
set.seed(123) # For reproducibility
model <- train(
  x = feature_data,
  y = days_to_death,
  method = "glmnet",
  tuneGrid = grid,
  trControl = control,
  metric = "RMSE" # Or another appropriate metric
)

# Display the best tuning parameters and model summary
print(model$bestTune)
print(model)

# Predict `days_to_death` using the best model
predictions <- predict(model, newdata = X_combined_new) # 'X_combined_new' is the new dataset for prediction

# Output predictions
print(predictions)

# Load necessary libraries
library(glmnet)
library(ggplot2)

# Assuming 'X_combined' is your feature matrix and 'days_to_death' is your response variable

# Fit the final model with the identified optimal parameters
final_model <- glmnet(feature_data, days_to_death, alpha = 0.25, lambda = 284.8036, family = "gaussian")

feature_data_matrix <- as.matrix(feature_data)

# Predict 'days_to_death' using the fitted model
predictions <- predict(final_model, newx = feature_data_matrix, s = 284.8036, type = "response")

# Convert predictions to a vector if it's not already
predictions <- as.vector(predictions)

# Create a data frame for plotting
plot_data <- data.frame(Actual = days_to_death, Predicted = predictions)

# Plot actual vs. predicted 'days_to_death'
ggplot(plot_data, aes(x = Actual, y = Predicted)) +
  geom_point(color = "blue") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  ggtitle("Actual vs. Predicted Days to Death") +
  xlab("Actual Days to Death") +
  ylab("Predicted Days to Death") +
  theme_minimal()

ggplot(data.frame(Days_to_Death = days_to_death), aes(x = Days_to_Death)) +
  geom_histogram(binwidth = 100, fill = "skyblue", color = "black") +
  ggtitle("Histogram of Actual Days to Death") +
  xlab("Days to Death") +
  ylab("Frequency") +
  theme_minimal()

# Combine actual and predicted data into a single data frame
data_for_plot <- data.frame(Actual = days_to_death, Predicted = predictions)
# Create the plot with actual vs. predicted values and a regression line
ggplot(data_for_plot, aes(x = Actual, y = Predicted)) +
  geom_point(aes(color = "Actual vs. Predicted"), alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "red") + # Adds a linear regression line
  labs(title = "Actual vs. Predicted Days to Death",
       x = "Actual Days to Death",
       y = "Predicted Days to Death",
       color = "Legend") +
  theme_minimal() +
  scale_color_manual(values = c("Actual vs. Predicted" = "blue", "Regression Line" = "red"))

# Using ggplot2 to plot actual vs. predicted values with a 1:1 line
library(ggplot2)

data_for_plot <- data.frame(Actual = days_to_death, Predicted = predictions)

ggplot(data_for_plot, aes(x = Actual, y = Predicted)) +
  geom_point(color = "blue", alpha = 0.5) +  # Plot actual vs. predicted values
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +  # Add a 1:1 line
  labs(title = "Actual vs. Predicted Days to Death",
       x = "Actual Days to Death",
       y = "Predicted Days to Death") +
  theme_minimal() +
  annotate("text", x = max(data_for_plot$Actual) * 0.8, y = max(data_for_plot$Predicted), label = "1:1 Line", color = "red")