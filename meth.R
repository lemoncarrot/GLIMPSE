library(GEOquery)
library(knitr)
library(limma)
library(minfi)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylationEPICmanifest)
library(RColorBrewer)
library(missMethyl)
library(minfiData)
library(Gviz)
library(DMRcate)
library(stringr)
library(GEOquery)
library(R.utils)
library(ggplot2)

df <- readRDS("GSE204943_meth.rds")
#cut out first row (IDs) ///////////////
df <- df[-1, ]
#second row is factors
#rename factors ///////////
df[1, ] <- sapply(df[1, ], function(x) {
  if(x == "sample_group: 1_Tumor") "tumor" else if(x == "sample_group: 2_OPL") "OPL" else if(x == "sample_group: 3_Normal") "normal" else x
})
cellTypes <- factor(df[1, ])
#Levels: sample_group: 1_Tumor sample_group: 2_OPL sample_group: 3_Normal
#additional filtering probably need to be performed here
df_t <- as.data.frame(t(df))
rownames(df_t) <- colnames(df)
colnames(df_t) <- rownames(df)
#df_t: samples as rownames, condition as first column, CpG sites rest of columns
samples_info <- data.frame(df_t[, 1])
colnames(samples_info) = "sample_type"
#column name "samples_info"
#rownames(samples_info) <- rownames(df_t)
#characteristics_ch1.1 - name of first column
#design needs to be rows as samples, conditions as columns (basically metadata)
design <- model.matrix(~0 + sample_type, data = samples_info)
#important! convert to numeric
df_t_numeric <- data.frame(lapply(df_t[, -1], function(x) as.numeric(as.character(x))))
rownames(df_t_numeric) <- rownames(df_t)
beta_values_matrix <- as.matrix(df_t_numeric[,-1])
beta_values_matrix <- t(beta_values_matrix)
#dimnames(beta_values_matrix) <- list(colnames(beta_values_matrix), rownames(beta_values_matrix))

#rows of design equals columns of beta_values_matrix
fit <- lmFit(beta_values_matrix, design)
contMatrix <- makeContrasts(sample_typetumor - sample_typenormal, levels=design)

fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)

# look at the numbers of DM CpGs at FDR < 0.05
summary(decideTests(fit2))

#adjust the column selection as per your actual annotation data
#ann850kSub <- ann850k[match(rownames(beta_values_matrix), ann850k$Name), 
#                      c("Name", "chr", "pos", "strand", "AddressA", "AddressB", "UCSC_RefGene_Name",
#                      "UCSC_RefGene_Accession", "UCSC_RefGene_Group")]
#ann850kSub <- ann850kSub[!is.na(rownames(ann850kSub)), ]

DMPs <- topTable(fit2, num=Inf, coef=1, genelist=ann850kSub)
head(DMPs)

dmr_site_names <- rownames(DMPs) #needs to be filtered
dmr_sites_df <- data.frame(SiteNames = dmr_site_names)
write.csv(dmr_sites_df, "dmr_site_names.csv", row.names = FALSE, quote = FALSE)
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
