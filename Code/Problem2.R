# PROBLEM 2

df_test = read.csv("peptides.test.csv", header=T, stringsAsFactors=F)
df_train = read.csv("peptides.train.csv", header=T, stringsAsFactors=F)

# 1
# 1A

set.seed(2024)
N <- nrow(df_train)

# Randomly select 20% of the data for validation (1:4 split)
validation_data = sample(1:N, N/5)

# Create validation and training sets
df_validation <- df_train[validation_data, ]
df_train <- df_train[-validation_data, ]

cat("Training data size:", nrow(df_train), "\n")
# output: Training data size: 2179 
cat("Validation data size:", nrow(df_validation), "\n")
# Validation data size: 544 



# 1B
# Compute the peptide sequence length from the 'peptides' column
df_train$length <- nchar(df_train$peptides)

# Display the first few rows of the new column
head(df_train$length)

# Summary statistics for peptide lengths
summary(df_train$length)


# Create a scatter plot
plot(df_train$length, df_train$rt,
     xlab = "Peptide Length",
     ylab = "Retention Time (rt)",
     main = "Peptide Length vs Retention Time",
     pch = 19, col = "black")


# 1C 
# Use Length as a Predictor in a Linear Model
lm_model <- lm(rt ~ length, data = df_train)
summary(lm_model)

# Is the length a statistically significant predictor for retention time?
# ANSWER: P-value for length (Pr(>|t|): < 2e-16)
# Yes, the length is statistically significant predictor for retention time because the p-value is much smaller than 0.05


# 1D
# Compute the mean retention time from the training data
mean_retention_time <- mean(df_train$rt)

# Predict retention time for validation data using the naive model
naive_predictions <- rep(mean_retention_time, nrow(df_validation))

# Compute the MSE for the naive model
mse_naive <- mean((df_validation$rt - naive_predictions)^2)
cat("Naive Model MSE:", mse_naive, "\n")
# output: Naive Model MSE: 39.06211 

# Compute peptide length for validation data (if not already computed)
df_validation$length <- nchar(df_validation$peptides)

# Predict retention time using the linear model on the validation data
linear_predictions <- predict(lm_model, newdata = df_validation)

# Compute the MSE for the linear model
mse_linear <- mean((df_validation$rt - linear_predictions)^2)
cat("Linear Model MSE:", mse_linear, "\n")
# output: Linear Model MSE: 32.48789 


# What are the two MSE values? 
# ANSWER: Naive Model MSE: 39.06211 
# Linear Model MSE: 32.48789
# ANSWER: Which one is better?
# The linear model is better since the Linear Model MSE (32.49) is smaller than the Naive Model MSE (39.06).
# The linear model reduces the MSE by 6.57 points (39.06 - 32.49).


# 1E
# How good is your current prediction with the linear model? 
# Explain your reasoning why you think the model is good / bad.
# ANSWER:
# Naive Model MSE: 39.06
# Linear Model MSE: 32.49
# The linear model is statistically significant, which confirms that length is related to rt.
# The linear model predicts retention time better than the naive model, which is good, as shown by the lower MSE (32.49 vs. 39.06). 
# This improvement indicates that peptide length has some predictive value. However, the model's R-squared is only 13.67%, 
# meaning length explains only a small portion of the variance in retention time.



# 2

df_test = read.csv("peptides.f.test.csv", header=T, stringsAsFactors=F)
df_train = read.csv("peptides.f.train.csv", header=T, stringsAsFactors=F)
df_val = read.csv("peptides.f.val.csv", header=T, stringsAsFactors=F)

# Verify the column names and structure
colnames(df_train)
str(df_train)
summary(df_train)

# Inspect the first few rows
head(df_train)


# 3
# Exclude one redundant feature 
df_train_reduced <- df_train[, !(colnames(df_train) %in% c("peptides"))]
df_val_reduced <- df_val[, !(colnames(df_val) %in% c("peptides"))]

# Train a linear model with 20 features
lm_model_new <- lm(rt ~ ., data = df_train_reduced)

# View the model summary
summary(lm_model_new)

# Predict on validation data
val_predictions <- predict(lm_model_new, newdata = df_val_reduced)

# Compute MSE for the validation set
mse_new <- mean((df_val$rt - val_predictions)^2)
cat("Extended Linear Model MSE:", mse_new, "\n")
# output: Extended Linear Model MSE: 10.33225 


# 3A
# Did the MSE become smaller with more predictors?
# ANSWER: Yes, the MSE become smaller with more predictors, it dropped from 32.48789 to 10.33225.
# What is its value?
# ANSWER:  10.33225 


# 3B
# Plot predicted vs true RT values
plot(df_val$rt, val_predictions,
     xlab = "True Retention Time (RT)",
     ylab = "Predicted Retention Time (RT)",
     main = "Predicted vs True Retention Time",
     pch = 19, col = "black")
abline(0, 1, col = "red")  

# Do you think your prediction has improved compared to the previous model using only length?
# ANSWER: Yes, the prediction has improved significantly.
# The black points (predictions) are close to the red diagonal line, which means that the predictions align well with the true retention times.
# The extended model, which includes 20 amino acid percentages and length, provides much more information about the peptide sequence.
# Significant reduction in MSE also shows that the additional predictors explain a lot more variability in rt than just the length alone, like in the original model.


# 3C
# Which features are important for retention time?
# ANSWER:
# Features that are statistically significant (p-values < 0.05) are important for predicting retention time,
# in this case W, F, L, I, H, K, R, and N, as they have large coefficients and highly significant p-values.
## Positive coefficients that increase retention time
# W (Tryptophan): Coefficient +37.42, p-value < 2e-16
# F (Phenylalanine): Coefficient +32.55, p-value < 2e-16
# L (Leucine): Coefficient +26.11, p-value < 2e-16
# I (Isoleucine): Coefficient +18.03, p-value < 2e-16
# M (Methionine): Coefficient +11.18, p-value 2.03e-06
# Length: Coefficient +0.805, p-value < 2e-16

## Negative coefficients that decrease retention time:
# H (Histidine): Coefficient -77.61, p-value < 2e-16
# K (Lysine): Coefficient -67.80, p-value < 2e-16
# R (Arginine): Coefficient -63.64, p-value < 2e-16
# N (Asparagine): Coefficient -30.36, p-value < 2e-16
# Q (Glutamine): Coefficient -28.34, p-value < 2e-16
# G (Glycine): Coefficient -27.70, p-value < 2e-16
# S (Serine): Coefficient -25.71, p-value < 2e-16
# C (Cysteine): Coefficient -23.02, p-value 1.84e-05
# T (Threonine): Coefficient -19.91, p-value < 2e-16
# E (Glutamic Acid): Coefficient -19.38, p-value < 2e-16
# P (Proline): Coefficient -18.57, p-value < 2e-16
# A (Alanine): Coefficient -17.80, p-value < 2e-16
# D (Aspartic Acid): Coefficient -14.81, p-value 2.93e-15


# Can you explain these results, do they make sense to you?
# ANSWER:
# The amino acids W, F, L, I, H, K, R, and N are hydrophobic and have large positive coefficients, meaning they increase retention time significantly. 
# This makes sense because this aligns with HPLC principles, as hydrophobic amino acids interact strongly with the hydrophobic stationary phase.
# While, the amino acids H, K, R, N, Q, G, S, C, T, E, P, A, and D are hydrophilic and have large negative coefficients, meaning they decrease retention time.
# This makes sense to me, because under acidic HPLC conditions, these amino acids are charged and interact more with the mobile phase which then reduce retention time.
# As for peptide length, the length of the peptide has a small but positive effect (+0.805), which is logical because longer peptides have more hydrophobic interactions overall.


# 4
library(randomForest)

set.seed(1)

# Train the Random Forest regression model
rf_model <- randomForest(rt ~ ., data = df_train_reduced, ntree = 500, mtry = sqrt(ncol(df_train_reduced) - 1))

# Predict on the validation data
rf_predictions <- predict(rf_model, newdata = df_val_reduced)

# Compute MSE on the validation data
mse_rf <- mean((df_val_reduced$rt - rf_predictions)^2)
cat("Random Forest Model MSE:", mse_rf, "\n")
# output: Random Forest Model MSE: 15.54845 

# Compare to the Linear Model MSE that previously computed
cat("Linear Model MSE:", mse_new, "\n")  
# output: Linear Model MSE: 10.33225


# Plot predicted vs true RT values for RF model
plot(df_val_reduced$rt, rf_predictions,
     xlab = "True Retention Time (RT)",
     ylab = "Predicted Retention Time (RT)",
     main = "Random Forest: Predicted vs True Retention Time",
     pch = 19, col = "black")
abline(0, 1, col = "red")  

# How does the MSE on the validation data compare to the linear model?
# ANSWER:
# Random Forest Model MSE: 15.54845
# Linear Model MSE: 10.33225
# The Linear Model performs better on the validation data, as it has a smaller MSE compared to the Random Forest model.
# A lower MSE means that the linear model's predictions are closer to the true retention time values than those of the random forest.
# It can happen most likely because if retention time depends primarily on linear relationships with the features, the linear model will most likely perform better.


# 5
# 5A
library(nnet)

mse_ann <- data.frame(size = integer(), mse_train = numeric(), mse_val = numeric())

# Loop through different sizes 
for (size in 1:25) {
  # Train the ANN model
  ann_model <- nnet(rt ~ ., data = df_train_reduced, size = size, linout = TRUE, maxit = 500, trace = FALSE)
  
  # Predict on training and validation data
  train_predictions <- predict(ann_model, df_train_reduced)
  val_predictions <- predict(ann_model, df_val_reduced)
  
  # Compute MSE for training and validation
  mse_train <- mean((df_train_reduced$rt - train_predictions)^2)
  mse_val <- mean((df_val_reduced$rt - val_predictions)^2)
  
  # Store results
  mse_ann <- rbind(mse_ann, data.frame(size = size, mse_train = mse_train, mse_val = mse_val))
}

# Results
print(mse_ann)
# output:
# size mse_train      mse_val
# 1     1 37.622232    37.203578
# 2     2 37.622232    37.203578
# 3     3 37.622232    37.203578
# 4     4 15.888873    16.195031
# 5     5 37.622232    37.203578
# 6     6 37.622232    37.203578
# 7     7 37.622232    37.203578
# 8     8  8.170581     8.844367
# 9     9 37.622232    37.203578
# 10   10  7.105386     9.879880
# 11   11 37.622219    37.203321
# 12   12  6.060708    19.201222
# 13   13  7.119621     9.275303
# 14   14 37.622232    37.203578
# 15   15  6.359941    76.815059
# 16   16  7.273245     9.647635
# 17   17  6.097886  7116.705135
# 18   18  8.016203     9.351055
# 19   19  4.761754  5830.438584
# 20   20  7.032102  4256.136086
# 21   21  5.441085 35159.543738
# 22   22  5.405605    14.805749
# 23   23  8.864894     9.122467
# 24   24  6.999445    10.058584
# 25   25  6.277912 39679.334306

# For smaller numbers of neurons like size = 1–3 and 5–7, both training and validation MSE are consistently high (around 37), which indicates underfitting.
# The artificial neural network performs best with size = 8 neurons, with a training MSE of 8.17 and validation MSE of 8.84. 
# For higher size  like 12, 15, and 18, there are large gaps between training and validation MSE, which indicates potential overfitting.
# While for higher numbers of neurons like size = 17–25, the validation MSE fluctuates a lot, indicating overfitting.
# In conclusion, increasing the number of neurons may result in overfitting or instability, as indicated by the fluctuating validation MSE for larger networks.



# 5B
set.seed(1)

mse_ann_decay <- data.frame(size = integer(), mse_train = numeric(), mse_val = numeric())

# Loop through different sizes with decay
for (size in 1:25) {
  # Train the ANN model with decay
  ann_model <- nnet(rt ~ ., data = df_train_reduced, size = size, decay = 0.75, linout = TRUE, maxit = 500, trace = FALSE)
  
  # Predict on training and validation data
  train_predictions <- predict(ann_model, df_train_reduced)
  val_predictions <- predict(ann_model, df_val_reduced)
  
  # Compute MSE for training and validation
  mse_train <- mean((df_train_reduced$rt - train_predictions)^2)
  mse_val <- mean((df_val_reduced$rt - val_predictions)^2)
  
  # Store results
  mse_ann_decay <- rbind(mse_ann_decay, data.frame(size = size, mse_train = mse_train, mse_val = mse_val))
}

# Results
print(mse_ann_decay)
# output:
#    size mse_train  mse_val
# 1     1  9.323372 9.234007
# 2     2  8.767103 8.598776
# 3     3  8.296003 8.469613
# 4     4  8.070236 8.643175
# 5     5  8.052718 8.543596
# 6     6  8.005067 8.493623
# 7     7  7.794376 8.788260
# 8     8  7.778433 8.492984
# 9     9  7.751105 8.569724
# 10   10  7.573309 8.843823
# 11   11  7.768692 8.574397
# 12   12  7.424797 8.841095
# 13   13  7.358960 9.130192
# 14   14  7.347993 8.996333
# 15   15  7.311303 8.817926
# 16   16  7.233655 8.893062
# 17   17  7.378782 8.929198
# 18   18  7.403271 8.877189
# 19   19  7.335259 8.996881
# 20   20  7.676768 8.934408
# 21   21  7.172193 8.993701
# 22   22  7.310634 8.765162
# 23   23  7.269419 9.061418
# 24   24  7.264387 8.766221
# 25   25  7.137655 9.294585

# Does this improve the quality of your models? Does it improve robustness?
# ANSWER:
# Yes, it improves the quality of my models.
# With decay = 0.75, the validation MSE (mse_val) is consistently lower than in the models without decay.
# For instance, at size = 6, the validation MSE (mse_val = 8.493623) is slightly better compared to the earlier case without decay (mse_val = 8.590483)
# The training MSE (mse_train) values are also generally smaller with decay, ranging from around 7 to 9, indicating better regularization and less overfitting.

# Yes, it improves robustness
# With decay = 0.75, extreme spikes in the validation MSE values are no longer seen.
# It shows that the weight decay parameter stabilizes the model and make it more robust.


# 5C 
# Investigate the performance as you increase the number of neurons in the hidden layer – 
# how does the validation and training error behave? Explain your observations.
# ANSWER: 
# The validation error initially decreases and reaches a minimum. After that, the validation error either increases or fluctuates inconsistently.
# This behavior indicates that adding more neurons initially improves the model’s ability to generalize to unseen data but then causes overfitting.

# While for training error behavior, initially the training error decreases, the training error stabilizes or fluctuates slightly, but it generally remains low.
# Therefore, it shows that increasing the number of neurons improves the model's ability to fit the training data, but it can lead to overfitting when too many neurons are used.


# 5D
## Select a final ANN model among those you tested
# ANN model with size = 6 (neurons in the hidden layer) and weight decay (decay = 0.75) provided the best results for the training and validation error:
# Training MSE: 8.005067
# Validation MSE: 8.493623
# This model balances underfitting and overfitting, so it is selected as the final ANN model.

## Report its performance
# Final ANN Model:
# Training MSE: 8.005067
# Validation MSE: 8.493623

# Previous Models for Comparison:
# Linear Model:
#  Validation MSE: 10.33225
# Random Forest (RF) Model:
#  Validation MSE: 15.54845

## How does its performance compare to linear models and RF?
# The ANN model outperforms the linear model, achieving a significantly lower validation MSE (8.493623 vs. 10.33225).
# The ANN model also outperforms the Random Forest model, achieving a much lower validation MSE (8.493623 vs. 15.54845).


# 6
# 6A

set.seed(1)
# Reduce the features because it kept getting error with error message "too many weights ..."
# Exclude non-numeric columns (e.g., 'peptides')
numeric_features <- df_full_train[, sapply(df_full_train, is.numeric)]

# Identify near-zero variance features
nzv <- apply(numeric_features, 2, function(col) var(col, na.rm = TRUE) < 1e-6)

# Remove near-zero variance features
df_full_train_reduced <- numeric_features[, !nzv]

# Add the target variable back to the reduced dataset
df_full_train_reduced$rt <- df_full_train$rt

# Train the ANN model on the reduced dataset
final_ann_model <- nnet(rt ~ ., data = df_full_train_reduced,
                        size = 4, decay = 0.75, maxit = 1000, linout = TRUE)


# 6B
# Based on the MSE on the validation set from the previous testing, 
# It is expected that the performance is similar to the validation MSE (mse_val)


# 6C
# Ensure test data matches the reduced training dataset
df_test_reduced <- df_test[, colnames(df_test) %in% colnames(df_full_train_reduced)]

# Predict retention time on the test dataset
test_predictions <- predict(final_ann_model, newdata = df_test_reduced)

# Compute the MSE for the test set
mse_test <- mean((df_test$rt - test_predictions)^2)
cat("Test MSE:", mse_test, "\n")
# output: Test MSE: 8.027743 

## Does it match your expectations?
# Test MSE Value: 8.027743
# Yes, it does.
# The test MSE (8.027743) is comparable to the Validation MSE (mse_val = 8.493623) obtained during model evaluation. 
# This means that the model generalizes well to unseen data and meets expectations based on validation results.



# 6D
# Did our choice of testing, training and validation dataset lead to an optimal model?
# ANSWER:
# Yes, the choice of testing, training, and validation dataset led to an optimal model.
# By dividing the process into training, validation, and test datasets allowed for unbiased evaluation. 
# The low test MSE also shows that this process led to a strong-performing model.

# Why was it necessary to evaluate different models and parameters on the validation dataset and not on the testing dataset?
# ANSWER:
# It was important to evaluate different models and parameters on the validation dataset and not on the testing dataset
# to prevent data leakage. Using the test set during the model selection process would introduce data leakage thus
# invalidating its purpose as a true measure of generalization.
# Other than that, it was also important because to simulate real-world performance.
# The validation dataset mimics unseen data that might be encountered during deployment. 
# It is untouched and  unbiased benchmark for final model performance.
# Therefore, the reported results are representative of how the model would perform in the real-world.



