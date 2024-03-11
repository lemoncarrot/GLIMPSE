library(R.utils)
library(limma)
library(data.table)
library(pheatmap)

#LOADING DATA
data <- fread("/home/dbr_journalclub/GSE90496_beta.txt.gz")
#cpg sites (450k) as rows and samples as columns
#pvalues included ("Detection Pval")
#class is "data.table" 
#cpg names in first column (not row names) ("ID_REF")
#beta values are already numeric
#column names are as follows: "SAMPLE X"

metadata <- read.csv("/home/dbr_journalclub/GSE90496_metadata.csv")
metadata <- metadata[, c("title", "geo_accession", "methylation.class.ch1")]
metadata$title <- sub(".*sample (\\d+).*", "SAMPLE \\1", metadata$title)
metadata <- as.data.table(metadata)
#dataframe with rows as samples and columns as characteristics
#"title" (SAMPLE X) "geo_accession" "methylation.class.ch1" ("GBM, ..." and "CONTR, ..." desired) 
#either filter now or later
cat("initial datasets loaded")

#DATA PREPROCESSING INTO FORMAT
pval_columns <- grep("Pval", names(data), value = FALSE)
data[, (pval_columns) := NULL]
metadata_transposed <- t(metadata[, -1, with = FALSE]) # Exclude the 'title' column for transposition
colnames(metadata_transposed) <- metadata$title
metadata_transposed_dt <- as.data.table(metadata_transposed)
setnames(metadata_transposed_dt, colnames(metadata_transposed_dt), paste0("SAMPLE ", 1:ncol(metadata_transposed_dt)))
metadata_for_merge <- data.table(ID_REF = c("geo_accession", "methylation.class.ch1"), metadata_transposed_dt)
cat("preprocessing done")
final_data <- rbindlist(list(metadata_for_merge, data), use.names = TRUE, fill = TRUE)
#dim(final_data)
#[1] 428801   2802
#class(final_data)
#[1] "data.table" "data.frame"
#no rownames
#colnames SAMPLE X
#first row GSM ID, second row glioma subtype identifier
#first column cpg site identifier

#DATA FILTERING FOR GBM AND CONTROL
final_data[2, ] <- lapply(final_data[2, ], function(x) {
  # Replace "ANYTEXT, XXX" with "ANYTEXT" (remove comma and anything after)
  gsub("([^,]*),.*", "\\1", x)
})

subtypes <- final_data[2, ]
selected_cols <- grepl("GBM|CONTR", subtypes)
df <- final_data[, c(TRUE, selected_cols), with=FALSE]

#dataframe of characters: converted to numeric late
#samples as columns and features as rows
#df <- readRDS("GSE204943_meth.rds")
#cut out first row (IDs) ///////////////
df <- df[-1, ]
#second row is factors
#rename factors ///////////
#next 3 lines is for dataframe
#df[1, ] <- sapply(df[1, ], function(x) {
#  if(x == "sample_group: 1_Tumor") "tumor" else if(x == "sample_group: 2_OPL") "OPL" else if(x == "sample_group: 3_Normal") "normal" else x
#})
#next 6 lines is for data.table
#might need to setDT
#run twice with control and normal
replace_function <- function(x) {
  x <- ifelse(grepl("GBM", x), "tumor", x)
  return(x)
}
df[1, (names(df)[-1]) := lapply(.SD, replace_function), .SDcols = names(df)[-1]]

#cellTypes <- factor(df[1, ])
#Levels: sample_group: 1_Tumor sample_group: 2_OPL sample_group: 3_Normal
#additional filtering probably need to be performed here
#df_t <- as.data.frame(t(df))
#rownames(df_t) <- colnames(df)
#colnames(df_t) <- rownames(df)
#df_t: samples as rownames, condition as first column, CpG sites rest of columns
#samples_info <- data.frame(df_t[, 1])
sample_names <- colnames(df)[-1]
cell_types <- df[1, -1]  
samples_info <- data.frame(CellType = t(cell_types), row.names = sample_names)
colnames(samples_info) = "sample_type"
#column name "samples_info"
#rownames(samples_info) <- rownames(df_t)
#characteristics_ch1.1 - name of first column
#design needs to be rows as samples, conditions as columns (basically metadata)
design <- model.matrix(~0 + sample_type, data = samples_info)
#important! convert to numeric
columns_to_modify <- setdiff(names(df), names(df)[1])
df[2:.N, (columns_to_modify) := lapply(.SD, function(x) as.numeric(as.character(x))), .SDcols = columns_to_modify]
beta_values_matrix <- as.matrix(df[-1, -1, with = FALSE])
rownames(beta_values_matrix) <- df$ID_REF[-1]
colnames(beta_values_matrix) <- colnames(df)[-1]
#now have numeric data table with first column as cpg sites and columns as samples 
#(filled in with beta values)
##df_t_numeric <- data.frame(lapply(df_t[, -1], function(x) as.numeric(as.character(x))))
##rownames(df_t_numeric) <- rownames(df_t)
##beta_values_matrix <- as.matrix(df_t_numeric[,-1])
##beta_values_matrix <- t(beta_values_matrix)

#dimnames(beta_values_matrix) <- list(colnames(beta_values_matrix), rownames(beta_values_matrix))

#rows of design equals columns of beta_values_matrix
numeric_mat <- matrix(as.numeric(as.character(beta_values_matrix)), nrow = nrow(beta_values_matrix), ncol = ncol(beta_values_matrix))
rownames(numeric_mat) <- rownames(beta_values_matrix)
colnames(numeric_mat) <- colnames(beta_values_matrix)
#saved numeric_mat rds
saveRDS(numeric_mat, "GSE90496_matrix_num2.rds")

fit <- lmFit(numeric_mat, design)
colnames(design) <- gsub(" ", "_", colnames(design)) #make syntactically valid names
colnames(design) <- gsub(",", "_", colnames(design)) #make syntactically valid names
colnames(design) <- gsub("/", "_", colnames(design)) #make syntactically valid names


contMatrix <- makeContrasts(sample_typetumor - sample_typenormal, levels=design)

fit2 <- contrasts.fit(fit, contMatrix)
#might be an issue
#Warning message:
#  In contrasts.fit(fit, contMatrix) :
#  row names of contrasts don't match col names of coefficients
#fit2 <- eBayes(fit2)

# look at the numbers of DM CpGs at FDR < 0.05
summary(decideTests(fit2))

#needs extra packages
#adjust the column selection as per your actual annotation data
#ann850kSub <- ann850k[match(rownames(beta_values_matrix), ann850k$Name), 
#                      c("Name", "chr", "pos", "strand", "AddressA", "AddressB", "UCSC_RefGene_Name",
#                      "UCSC_RefGene_Accession", "UCSC_RefGene_Group")]
#ann850kSub <- ann850kSub[!is.na(rownames(ann850kSub)), ]

DMPs <- topTable(fit2, num=Inf, coef=1)
head(DMPs)
filtered_DMPs <- DMPs[abs(DMPs$logFC) > 0.1 & DMPs$adj.P.Val < 0.05, ]
head(filtered_DMPs)
dim(filtered_DMPs)

dmr_site_names <- rownames(filtered_DMPs) #needs to be filtered
dmr_sites_df <- data.frame(SiteNames = dmr_site_names)
write.csv(dmr_sites_df, "DMR_0_1.csv", row.names = FALSE, quote = FALSE)
#KEGG analysis after this

#volcano
#can adjust filtering 
png(filename="meth_1.png", width=8, height=6, units="in", res=300)
ggplot(DMPs, aes(x = logFC, y = -log10(P.Value), color = (adj.P.Val < 0.05 & abs(logFC)>0.1))) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("grey", "red")) +
  labs(title = "Volcano Plot of DMPs",
       x = "Log Fold Change",
       y = "-Log10 P-Value") +
  theme_minimal() +
  theme(legend.title=element_blank())
dev.off()

png(filename="meth_2.png", width=8, height=6, units="in", res=300)
ggplot(DMPs, aes(x = logFC)) +
  geom_histogram(bins = 50, fill = "red", color = "white") +
  labs(title = "Distribution of Log Fold Changes",
       x = "Log Fold Change",
       y = "Frequency") +
  theme_minimal() +
  theme(legend.title=element_blank())
dev.off()

#heatmap here
#possible circular graph here
#needs to integrate ann850k with beta values, p-values, and sample conditions
#needs adjusting

sorted_DMPs <- filtered_DMPs[order(-abs(filtered_DMPs$logFC)), ]
top100_DMPs <- head(sorted_DMPs, 100)
top100_sites <- rownames(top100_DMPs)
cpg_annotations <- data.frame(
  logFC = DMPs$logFC,
  adj.P.Value = DMPs$adj.P.Val
)
rownames(cpg_annotations) <- rownames(numeric_mat)
#samples_info = sample annotation
dmr_mat <- numeric_mat[top100_sites, ]
cpg_annotations <- cpg_annotations[top100_sites, ]

#cool one
png(filename="cool_3.png", width=16, height=10, units="in", res=300)
pheatmap(dmr_mat,
         annotation_col = samples_info,    # Column annotations: Sample attributes like tumor/normal
         annotation_row = cpg_annotations,      # Row annotations: CpG site attributes like p-values
         scale = "row",                        
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = TRUE,
         show_colnames = FALSE)
dev.off()

tumor_normal_indices <- which(samples_info$sample_type %in% c("tumor", "normal"))
samplefilter <- rownames(samples_info)[c(tumor_normal_indices)]
dmr_mat2 <- numeric_mat[top100_sites, samplefilter]

#real one
png(filename="meth_3.png", width=16, height=10, units="in", res=300)
pheatmap(dmr_mat2,
         annotation_col = samples_info,    # Column annotations: Sample attributes like tumor/normal
         annotation_row = cpg_annotations,      # Row annotations: CpG site attributes like p-values
         scale = "row",                       
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = TRUE,
         show_colnames = FALSE)
dev.off()

