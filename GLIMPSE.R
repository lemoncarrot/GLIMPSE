library(reshape2)
library(readr) 


# Step 1: Read the CSV file containing the site names
cg_sites_df <- read.csv("DMR_filtered_overlap.csv")
cg_site_names <- cg_sites_df$x

# Step 2: Load the TCGA-GBM RDS file into a dataframe
tcga_gbm_df <- readRDS("TCGA-GBM_meth.rds")

# Step 3: Ensure that your filtered dataset includes the cg sites as well as 'sample' and 'days_to_death'
# Prepare a vector of columns to keep, making sure to include your metadata columns
columns_to_keep <- c(cg_site_names, "sample", "days_to_death", "vital_status")

# Filter the dataset to include only the specified columns, ensuring to check for column existence
filtered_tcga_gbm_df <- tcga_gbm_df[, colnames(tcga_gbm_df) %in% columns_to_keep]
filtered_tcga_gbm_df <- filtered_tcga_gbm_df[, colSums(is.na(filtered_tcga_gbm_df)) < nrow(filtered_tcga_gbm_df)]
filtered_tcga_gbm_df <- filtered_tcga_gbm_df[!is.na(filtered_tcga_gbm_df$days_to_death), ]


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
library(survminer) # For ggsurvplot
library(timeROC)
library(glmnetUtils)
library(FDb.InfiniumMethylation.hg19)
library(randomForest)
library(boot)



#-------
#LOAD DATA
sig_genes <- read.csv('DMR_DEG_overlap_genes.csv')
tcga_rna <- readRDS('TCGA-GBM_rna.rds')
sig_genes <- unlist(sig_genes)

library(dplyr)
library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genes <- getBM(filters="hgnc_symbol", 
               attributes=c("ensembl_gene_id","hgnc_symbol"), 
               values=sig_genes, 
               mart=ensembl)
mapping_vector <- setNames(genes$ensembl_gene_id, genes$hgnc_symbol)
ensembl_ids <- mapping_vector[sig_genes]
ensembl_ids <- as.data.frame(ensembl_ids)
#FORMAT DATA
# Convert to list format and filter data
gene_ids <- ensembl_ids$ensembl_ids
colnames(tcga_rna) <- sapply(strsplit(colnames(tcga_rna), "\\."), `[`, 1)
filtered_rna_data <- tcga_rna[, c("sample", "days_to_death", "vital_status", grep(paste0("^", paste(gene_ids, collapse="|"), "$"), colnames(tcga_rna), value = TRUE))]

# Filter for samples overlapping in RNA and methylation data
overlap_samples <- intersect(filtered_rna_data$sample, filtered_tcga_gbm_df$sample)
filtered_rna_data <- filtered_rna_data[filtered_rna_data$sample %in% overlap_samples, ]
filtered_tcga_gbm_df <- filtered_tcga_gbm_df[filtered_tcga_gbm_df$sample %in% overlap_samples, ]

# Ensure no missing values in survival data and combine RNA and methylation data
complete_cases <- complete.cases(filtered_rna_data$days_to_death, filtered_rna_data$vital_status)
filtered_rna_data <- filtered_rna_data[complete_cases, ]
filtered_tcga_gbm_df <- filtered_tcga_gbm_df[complete_cases, ]
filtered_tcga_gbm_df_no_na <- na.omit(filtered_tcga_gbm_df)

combined_data <- cbind(filtered_rna_data, filtered_tcga_gbm_df[,-c(1:3)]) # Avoid duplicating 'sample', 'days_to_death', and 'vital_status' columns
combined_data_filtered <- combined_data[, colSums(is.na(combined_data)) < nrow(combined_data)]
complete_cols <- colSums(is.na(combined_data_filtered)) == 0
data_complete_cols <- combined_data_filtered[, complete_cols]
setDT(data_complete_cols)

# Ensure the samples overlap after filtering by days_to_death
overlap_samples <- intersect(filtered_rna_data$sample, filtered_tcga_gbm_df$sample)
filtered_rna_data <- filtered_rna_data[filtered_rna_data$sample %in% overlap_samples, ]
filtered_tcga_gbm_df <- filtered_tcga_gbm_df[filtered_tcga_gbm_df$sample %in% overlap_samples, ]

# Continue with creating the combined dataset and proceeding with your analysis
combined_data <- cbind(filtered_rna_data, filtered_tcga_gbm_df[,-c(1:3)]) # Avoid duplicating 'sample', 'days_to_death', and 'vital_status' columns
combined_data_filtered <- combined_data[, colSums(is.na(combined_data)) < nrow(combined_data)]
complete_cols <- colSums(is.na(combined_data_filtered)) == 0
data_complete_cols <- combined_data_filtered[, complete_cols]
setDT(data_complete_cols)


# At this point, both filtered_rna_data and filtered_tcga_gbm_df have been filtered to exclude samples with days_to_death over your cutoff
# You can now proceed with your model training and testing as planned


columns_of_interest <- grep("^cg|sample|days_to_death|vital_status$", colnames(filtered_tcga_gbm_df), value = TRUE)

# Subset the DataFrame to keep only the selected columns
filtered_tcga_gbm_df <- filtered_tcga_gbm_df[, columns_of_interest]


# Extract survival-related data into a separate dataframe
survival_data <- data_complete_cols[, .(days_to_death, vital_status)]
feature_data <- data_complete_cols[, !names(data_complete_cols) %in% c("sample", "days_to_death", "vital_status"), with = FALSE]
days_to_death <- data_complete_cols$days_to_death
vital_status_binary <- as.integer(survival_data$vital_status == "Dead")
survival_data$vital_status_binary <- as.integer(survival_data$vital_status == "Dead")


# Convert data.table to matrix if necessary
if (is.data.table(feature_data)) {
  feature_data_matrix <- as.matrix(feature_data)
} else {
  feature_data_matrix <- feature_data # Assuming feature_data is already a matrix
}

# PCA with variance explained criterion
preProcVar <- preProcess(feature_data_matrix, method = "pca", thresh = 0.9) # Adjust thresh for desired variance explained
feature_data_pca <- predict(preProcVar, feature_data_matrix)
library(caret)

# Update the feature dataset
feature_data <- feature_data_pca

# Convert the matrix to a data frame for easier manipulation
feature_data_df <- as.data.frame(feature_data_matrix)

# Initialize a list to hold all the polynomial features
poly_features <- list()

# Generate polynomial features for each feature
for(feature_name in names(feature_data_df)) {
  # Second degree polynomial features
  poly_features[[paste(feature_name, "squared", sep = "_")]] <- feature_data_df[[feature_name]]^2
  
  # Interaction terms can also be added here if needed
}

# Convert the list of polynomial features into a dataframe
poly_features_df <- as.data.frame(poly_features)

# Bind the original features with the polynomial features
feature_data_enhanced <- cbind(feature_data_df, poly_features_df)

# If needed, convert back to a matrix for modeling
feature_data <- as.matrix(feature_data_enhanced)

# Create the survival object
surv_obj <- Surv(time = survival_data$days_to_death, event = vital_status_binary)


# Extract time and status from surv_obj (assuming surv_obj can be treated like a dataframe for extraction, which it cannot directly. This is for conceptual understanding)
time <- survival_data$days_to_death  # Time to event or censoring
status <- vital_status_binary  # Event occurred (1) or censored (0)
y <- matrix(c(time, status), ncol = 2, byrow = FALSE)
colnames(y) <- c("time", "status")

x <- as.matrix(feature_data)

#cv and gridsearch
control <- trainControl(method = "cv", number = 5, search = "grid")
grid <- expand.grid(
  alpha = seq(0, 1, length.out = 5), # Adjust the sequence as needed
  lambda = 10^seq(-3, 3, length.out = 100) # Adjust the sequence as needed
)
# Refining the grid for hyperparameter tuning
refined_grid <- expand.grid(alpha = seq(0, 1, by = 0.1), # More refined alpha grid
                            lambda = 10^seq(-3, 1, by = 0.2)) # More refined lambda grid
set.seed(123) # For reproducibility
model <- train(
  x = feature_data,
  y = days_to_death,
  method = "glmnet",
  tuneGrid = refined_grid,
  trControl = control,
  metric = "RMSE" # Or another appropriate metric
)

# Display the best tuning parameters and model summary
print(model$bestTune)
print(model)


# Automatically extract the optimal alpha and lambda values
optimal_alpha <- model$bestTune$alpha
optimal_lambda <- model$bestTune$lambda

# Convert feature_data to a matrix, as required by glmnet
feature_data_matrix <- as.matrix(feature_data)

# Fit the final model with the identified optimal parameters
final_model <- glmnet(feature_data_matrix, days_to_death, alpha = optimal_alpha, lambda = optimal_lambda, family = "gaussian")

# Predict 'days_to_death' using the fitted model with the optimal lambda
# Note: The 's' parameter in the predict function is used to specify the value of lambda for prediction
predictions <- predict(final_model, newx = feature_data_matrix, s = optimal_lambda, type = "response")

# Convert predictions to a vector if it's not already
predictions <- as.vector(predictions)

# Create a data frame for plotting
plot_data <- data.frame(Actual = days_to_death, Predicted = predictions)


#-----
library(randomForest)

# Prepare the dataset (assuming 'feature_data' and 'days_to_death' are ready)
set.seed(123)
rf_model <- randomForest(x = feature_data_matrix, y = days_to_death)

# Predict using the Random Forest model
rf_predictions <- predict(rf_model, feature_data_matrix)

# Calculate and print the correlation coefficient
rf_cor_coef <- cor(rf_predictions, days_to_death)
print(paste("Random Forest Correlation Coefficient:", rf_cor_coef))

# Calculating RMSE
rmse_rf <- sqrt(mean((rf_predictions - days_to_death)^2))
print(paste("Random Forest RMSE:", rmse_rf))

# Plotting Actual vs. Predicted values
plot_data_rf <- data.frame(Actual = days_to_death, Predicted = rf_predictions)
ggplot(plot_data_rf, aes(x = Actual, y = Predicted)) +
  geom_point(alpha = 0.5) +
  geom_abline(color = "red", linetype = "dashed") +
  labs(title = "Random Forest: Actual vs. Predicted Days to Death", x = "Actual Days to Death", y = "Predicted Days to Death") +
  theme_minimal()

#-----
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


# Calculate the correlation coefficient
correlation_coefficient <- cor(days_to_death, predictions)

# Calculate RMSE
rmse <- sqrt(mean((days_to_death - predictions)^2))

# Print the results
print(paste("Correlation Coefficient:", correlation_coefficient))
print(paste("RMSE (in days):", rmse))


#-----------------------------


