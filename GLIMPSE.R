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
sig_genes <- read.csv('temp_67_genes.csv')
tcga_rna <- readRDS('TCGA-GBM_rna.rds')
tcga_meth <- readRDS('TCGA-GBM_meth.rds')

# Extract significant gene IDs
gene_ids <- sig_genes$Gene
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


# Calculate the correlation coefficient
correlation_coefficient <- cor(days_to_death, predictions)

# Calculate RMSE
rmse <- sqrt(mean((days_to_death - predictions)^2))

# Print the results
print(paste("Correlation Coefficient:", correlation_coefficient))
print(paste("RMSE (in days):", rmse))
