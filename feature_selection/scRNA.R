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

x <- readRDS("GSE84465.rds") #scrna dataset
gene_data <- as.data.frame(x[60:nrow(x), ])
metadata <- x[1:59, ]
metadata <- t(metadata)
metadata <- as.data.frame(metadata)
seurat_obj <- CreateSeuratObject(counts = gene_data, project = "GBM-scRNA", meta.data = metadata)

tumor_obj <- subset(seurat_obj, subset = characteristics_ch1.3 == "tissue: Tumor")
peripheral_obj <- subset(seurat_obj, subset = characteristics_ch1.3 == "tissue: Periphery")

#visualizing example DEGs in scRNA
seurat_obj <- AddMetaData(seurat_obj, metadata = seurat_obj@meta.data$characteristics_ch1.3, col.name = "tissue_type")
Idents(seurat_obj) <- seurat_obj@meta.data$characteristics_ch1.3
degs <- FindMarkers(seurat_obj, ident.1 = "tissue: Tumor", ident.2 = "tissue: Periphery", logfc.threshold = 0.25)
seurat_obj$active.ident <- seurat_obj$tissue_type
seurat_obj <- subset(seurat_obj, min.cells = 3, min.features=100)
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 1500)
seurat_obj <- ScaleData(seurat_obj)
degs <- FindMarkers(seurat_obj, ident.1 = "tissue: Tumor", ident.2 = "tissue: Periphery", logfc.threshold = 0.25)
top20_degs <- head(degs[order(degs$p_val_adj), ], 50)
top20_genes <- rownames(top20_degs)
png(filename="scRNA_9.png", width=22, height=10, units="in", res=300)
DoHeatmap(seurat_obj, features = top20_genes) + NoLegend()
dev.off()
degs$significant <- degs$p_val_adj < 0.05 & abs(degs$avg_log2FC) > 2
png(filename="scRNA_10.png", width=8, height=6, units="in", res=300)
ggplot(degs, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significant)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("black", "red")) +
  labs(title = " ",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") +
  theme_minimal() +
  theme(legend.position = "none")
dev.off()

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
#pdf("scRNA_4.pdf", width = 16, height = 10)
#DoHeatmap(tumor_obj, features = marker_genes)
#dev.off()
unique_genes <- unique(significant_tumor_markers$gene)
genes_df <- data.frame(Gene = unique_genes)
write.csv(genes_df, "scRNA_t_genes.csv", row.names = FALSE)

# Process Peripheral Data
peripheral_obj <- subset(peripheral_obj, min.cells = 3, min.features=100)
peripheral_obj[["percent.mt"]] <- PercentageFeatureSet(peripheral_obj, pattern = "^MT-")
peripheral_obj <- subset(peripheral_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)
#ncol(peripheral_obj)
#pdf("scRNA_5.pdf", width = 8, height = 6)
#VlnPlot(peripheral_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#dev.off()
peripheral_obj <- NormalizeData(peripheral_obj, normalization.method = "LogNormalize", scale.factor = 10000)
peripheral_obj <- FindVariableFeatures(peripheral_obj, selection.method = "vst", nfeatures = 1500)
#pdf("scRNA_6.pdf", width = 8, height = 6)
#VariableFeaturePlot(peripheral_obj)
#dev.off()
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
png(filename = "scRNA_7.png", width = 8, height = 6, units="in", res=300)
DimPlot(peripheral_obj, reduction = "tsne", group.by = "SingleR.celltype_adjusted", label = TRUE)
dev.off()
peripheral_markers <- FindAllMarkers(peripheral_obj, min.pct = 0.25, logfc.threshold = 0.25)
significant_peripheral_markers <- subset(peripheral_markers, p_val_adj < 0.05)
top10_markers2 <- significant_peripheral_markers %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)
marker_genes2 <- top10_markers2$gene
#pdf("scRNA_8.pdf", width = 16, height = 10)
#DoHeatmap(peripheral_obj, features = marker_genes)
#dev.off()
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
top_genes_table <- as.data.frame(top_genes)
gene_names <- rownames(top_genes_table)
gene_names_df <- data.frame(Gene = gene_names)
write.csv(gene_names_df, "scRNA_both_genes.csv", row.names = FALSE, quote = FALSE)
#2269 genes
saveRDS(lrt, file="lrt_results.rds")
