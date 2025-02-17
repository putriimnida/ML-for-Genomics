# PROBLEM 1

set.seed(123)
df = read.csv("raw_data_accuracy_benchmark.fix.txt", sep="\t")

N = length(unique(df$Sequence))
test_seq = sample(unique(df$Sequence), N/5)
df_test = subset(df, Sequence %in% test_seq)
df_train = subset(df, !(Sequence %in% test_seq))


# 1

df$log_oligo_length <- log(df$X.oligo.)
df_test$log_oligo_length <- log(df_test$X.oligo.)
df_train$log_oligo_length <- log(df_train$X.oligo.)

# linear model
model <- lm(Tm_exp ~ X.CG + Length + X.salt. + X.oligo. + log_oligo_length, data = df_train)

pred_training <- predict(model, df_train)
pred_testing <- predict(model, df_test)

mse_training <- mean((df_train$Tm_exp - pred_training)^2)
mse_testing <- mean((df_test$Tm_exp - pred_testing)^2)

cat("Training MSE:", mse_training, "\n") 
# output: Training MSE: 8.365003 
cat("Testing MSE:", mse_testing, "\n")
# output: Testing MSE: 7.027037 

summary(model)

# QUESTION: Based on the linear model, which variables are important for the prediction?
# ANSWER: 
# The most important predictors of Tm_exp are salt concentration, length, and GC content, 
# as shown by their small p-values and significant coefficients. The other two which are 
# Log_oligo_length is moderately significant, while X.oligo (oligo concentration) is not statistically significant.
# More explanation:
## X.CG (GC content): a large t-value (38.407) and a very small p-value (2e-16), 
# indicating that GC content is a critical predictor of Tm_exp.
## Length: 1 unit increase in sequence length increases Tm_exp by 0.8837. 
# Its significance (19.509 t-value) and very small p-value (< 2e-16)
## X.salt (salt concentration): the most important, a very large coefficient (12.02) and very small p-value (< 2e-16),
# meaning that salt concentration has a strong effect on the melting temperature.


# 2 
# create scatterplot
plot(df_train$Length, df_train$log_oligo_length,
     xlab = "Length",
     ylab = "Logarithm of Length",
     main = "Scatter Plot: Length vs. Logarithm of Length",
     pch = 20, col = "black")

# Interpretation:
## What do you see?
# The scatter plot shows the relationship between Length (x-axis) and Logarithm of Length (y-axis). 
# Each point represents a data entry which appear to be a clear, non-random pattern, suggesting some level of relationship between the two variables.

## What kind of relationship/model?
# There seems to be a linear relationship between the length and its logarithmic transformation.

## Are they correlated?
# The correlation coefficient between Length and Logarithm of Length is -0.242, which indicates a weak negative correlation. 
# This means that as Length increases, Logarithm of Length tends to decrease slightly, yet the relationship is not strong.
cor(df_train$Length, df_train$log_oligo_length)
# output: -0.2420364


# model A (contains the length but not the logarithm of the length)
model_a <- lm(Tm_exp ~ X.CG + Length + X.salt. + X.oligo., data = df_train)
pred_test_a <- predict(model_a, df_test)
mse_test_a <- mean((df_test$Tm_exp - pred_test_a)^2)
cat("MSE for model with Length only:", mse_test_a, "\n")
# output: MSE for model with Length only: 7.10059 

# model B (contains the logarithm of the length but not the length)
model_b <- lm(Tm_exp ~ X.CG + log_oligo_length + X.salt. + X.oligo., data = df_train)
pred_test_b <- predict(model_b, df_test)
mse_test_b <- mean((df_test$Tm_exp - pred_test_b)^2)
cat("MSE for model with Logarithm of Length only:", mse_test_b, "\n")
# output: MSE for model with Logarithm of Length only: 20.89625 

# model C (contains both)
model_c <- lm(Tm_exp ~ X.CG + Length + log_oligo_length + X.salt. + X.oligo., data = df_train)
pred_test_c <- predict(model_c, df_test)
mse_test_c <- mean((df_test$Tm_exp - pred_test_c)^2)
cat("MSE for model with both Length and Logarithm of Length:", mse_test_c, "\n")
# output: MSE for model with both Length and Logarithm of Length: 7.027037 

# Which performs the best?
# Model C (with both Length and Logarithm of Length) performs the best because it achieves the lowest MSE, 
# which means it makes the most accurate predictions on the test set.


# 3
install.packages("randomForest")
library(randomForest)

# a) train RF on the three models
# RF on model A
rf_a <- randomForest(Tm_exp ~ X.CG + Length + X.salt. + X.oligo., data = df_train)
pred_test_rf_a <- predict(rf_a, df_test)
mse_rf_a <- mean((df_test$Tm_exp - pred_test_rf_a)^2)
cat("RF MSE for model with Length only:", mse_rf_a, "\n")

# RF on model B
rf_b <- randomForest(Tm_exp ~ X.CG + log_oligo_length + X.salt. + X.oligo., data = df_train)
pred_test_rf_b <- predict(rf_b, df_test)
mse_rf_b <- mean((df_test$Tm_exp - pred_test_rf_b)^2)
cat("RF MSE for model with Logarithm of Length only:", mse_rf_b, "\n")

# RF on model C
rf_c <- randomForest(Tm_exp ~ X.CG + Length + log_oligo_length + X.salt. + X.oligo., data = df_train)
pred_test_rf_c <- predict(rf_c, df_test)
mse_rf_c <- mean((df_test$Tm_exp - pred_test_rf_c)^2)
cat("RF MSE for model with both Length and Logarithm of Length:", mse_rf_c, "\n")

# b) do you see the same results? Which models work best?
# ANSWER: Yes, I see the same result with the linear models. Model C achieves the lowest MSE (7.027037), making it the best model.

# Create variable of importance plot
varImpPlot(rf_a, main = "Variable Importance: RF Model (Length Only)")
varImpPlot(rf_b, main = "Variable Importance: RF Model (Log Length Only)")
varImpPlot(rf_c, main = "Variable Importance: RF Model (Both Length and Log Length)")

# c) Now compare the variable importance plots (varImpPlot) of the two models with four variables, what do you notice?
# comparsion:
# Variable Importance for RF Model with Length Only (model A):
# X.CG and X.salt. appear to have high importance (IncNodePurity).
# Length also contributes meaningfully but is slightly less important than X.CG and X.salt..
# X.oligo. shows minimal importance.

# Variable Importance for RF Model with Logarithm of Length Only (model B):
# Similar to the model A, X.CG and X.salt. remain the most important variables.
# log_oligo_length is less important compared to X.CG and X.salt., similar to how Length performed in the first model.
# X.oligo. also shows minimal importance.

# Conclusion:
# Regardless which one to remove, either length or log_oligo_length, X.CG and X.salt. are consistently the most important variables in both models.


# (d) Now compare the variable importance plots with the output of the linear model. 
# Which variables are important for the linear model, which ones are important for the RF?
# X.CG is highly significant with a strong effect on melting temperature.
# X.salt.is also highly significant, with a large coefficient indicating strong influence.
# Length is significant, though its contribution is smaller compared to X.CG and X.salt.
# log_oligo_length is moderately significant in the combined model but not as strong as Length.
# Both X.CG and X.salt are also important for the RF.
# However, when both length and oligo_length are included in Model C, RF assigns higher importance to Length compared to log_oligo_length.
# Thus, this aligns with the linear model’s results, where Length was more significant than log_oligo_length.
# Are the same variables important?
# Both X.CG and X.salt. are the most important predictors for Tm_exp in linear models and RF.


# 4
install.packages("nnet")
library(nnet)

# 4A
set.seed(3)

# Model A (Length only, size = 5)
ann_a <- nnet(Tm_exp ~ X.CG + Length + X.salt. + X.oligo., 
              data = df_train, 
              size = 5, decay = 0.1, maxit = 1000, linout = TRUE)
pred_train_a <- predict(ann_a, df_train)
pred_test_a <- predict(ann_a, df_test)
mse_train_a <- mean((df_train$Tm_exp - pred_train_a)^2)
mse_test_a <- mean((df_test$Tm_exp - pred_test_a)^2)

# Model B (log_oligo_length only, size = 5)
ann_b <- nnet(Tm_exp ~ X.CG + log_oligo_length + X.salt. + X.oligo., 
              data = df_train, 
              size = 5, decay = 0.1, maxit = 1000, linout = TRUE)
pred_train_b <- predict(ann_b, df_train)
pred_test_b <- predict(ann_b, df_test)
mse_train_b <- mean((df_train$Tm_exp - pred_train_b)^2)
mse_test_b <- mean((df_test$Tm_exp - pred_test_b)^2)

# Model C (Both Length and log_oligo_length, size = 5)
ann_c <- nnet(Tm_exp ~ X.CG + Length + log_oligo_length + X.salt. + X.oligo., 
              data = df_train, 
              size = 5, decay = 0.1, maxit = 1000, linout = TRUE)
pred_train_c <- predict(ann_c, df_train)
pred_test_c <- predict(ann_c, df_test)
mse_train_c <- mean((df_train$Tm_exp - pred_train_c)^2)
mse_test_c <- mean((df_test$Tm_exp - pred_test_c)^2)

# Print MSE results
cat("Model A - Train MSE:", mse_train_a, "Test MSE:", mse_test_a, "\n")
# output seed(1): Model A - Train MSE: 3.381881 Test MSE: 4.265583 
# output seed(2): Model A - Train MSE: 5.229487 Test MSE: 4.620813 
# output seed(3): Model A - Train MSE: 3.372627 Test MSE: 2.719141 
cat("Model B - Train MSE:", mse_train_b, "Test MSE:", mse_test_b, "\n")
# output seed(1): Model B - Train MSE: 13.89047 Test MSE: 15.11639 
# output seed(2): Model B - Train MSE: 14.45938 Test MSE: 15.05877
# output seed(3): Model B - Train MSE: 12.54343 Test MSE: 15.72506 
cat("Model C - Train MSE:", mse_train_c, "Test MSE:", mse_test_c, "\n")
# output seed(1): Model C - Train MSE: 4.218148 Test MSE: 4.534806 
# output seed(2): Model C - Train MSE: 2.762754 Test MSE: 2.296193
# output seed(3): Model C - Train MSE: 2.670827 Test MSE: 2.295687 


## What do you notice when you train the same three models? Do you see a difference? 
# Model C achieves the lowest train and test MSE, which shows that both Length and log_oligo_length improves prediction accuracy.
# Model B has higher train and test MSE, showing that log_oligo_length alone is insufficient to explain the variability in Tm_exp.
# Model A performs better than Model B but not as well as Model C

## When you repeat the training, does the model performance vary?
# Yes, when I repeat the training process with three different seeds, it leads to minor changes in MSE values.
# However, Model C remains the model with the lowest train and test MSE, followed by Model A, and Model C.

# 4B
set.seed(123)
# Train ANN with new parameters
ann_large <- nnet(Tm_exp ~ X.CG + Length + log_oligo_length + X.salt. + X.oligo., data = df_train,
                  size = 50, decay = 0.5, maxit = 5000, linout = TRUE)

pred_train_large <- predict(ann_large, df_train)
pred_test_large <- predict(ann_large, df_test)

# Compute the MSE
mse_train_large <- mean((df_train$Tm_exp - pred_train_large)^2)
mse_test_large <- mean((df_test$Tm_exp - pred_test_large)^2)

# Print 
cat("MSE for Large ANN - Train:", mse_train_large, "\n")
# output: MSE for Large ANN - Train: 1.545923 
cat("MSE for Large ANN - Test:", mse_test_large, "\n")
# output: MSE for Large ANN - Test: 2.795042 

# What do you observe? Do you see a difference in the MSE, do you see a difference between the test and training data MSE – why do you think that is?
# MSE for Large ANN - Train: 1.545923 
# MSE for Large ANN - Test: 2.795042 
# The relatively close values (1.545923 vs. 2.795042) suggest that the model generalizes reasonably well for a larger ANN.
# The larger ANN has a slightly better Test MSE compared to simpler ANN models in 4A

# 4C
# Based on your results, would you recommend ANN over RF or linear models for this problem?
# Is it a good idea to evaluate multiple models on the test data and select the best model based on its performance on the test set? 
# Explain your reasoning.

# ANSWER:
# We can use RF if the dataset size is small, and interpretability is important.
# But we can use ANN if we want to achieve the lowest MS and if we can ensure regularization to avoid overfitting.

# explanation:
## ANN has strengths:
# 1. Ability to Model Nonlinear Relationships:
# ANNs, especially larger ones, can capture complex relationships in the data that may not be easily captured by linear models or Random Forest (RF).
# In this case, the large ANN (size = 50) performed well, achieving a test MSE of 2.795042, close to the best simpler ANN model in 4A (2.498492).
# 2. Better Fit to Data:
# The ANN achieved the lowest train MSE (~1.5) compared to RF and linear models. This indicates its ability to fit the data well.
# ANN has limitations:
# 1. Risk of overfitting
# The initial results (before setting a seed) showed overfitting, with a high test MSE (4.071019).
# 2. The size of data
# ANNs often require large datasets to fully leverage their capabilities. With smaller datasets, simpler models like RF or linear regression are often more robust.
## RF Pros:
# RF achieved competitive test MSE (likely close to ~2.5 based on earlier results) and is less sensitive to overfitting.
## Linear model pros:
# The linear model with both Length and log_oligo_length achieved competitive MSE (7.027037), though not as good as RF or ANN.

# ANSWER:
# No, it is not a good idea. 
# This is because the model was fitted to the training data and the parameters were optimized to reduce the error on the training data. 
# The testing data represents 'new data,' which the model has not seen before, and is used to evaluate how well the model generalizes to unseen data. 
# By evaluating multiple models on the test set and selecting the best one, we indirectly risk tuning the model to the test set, 
# making it no longer an unbiased evaluation of performance on unseen data. 
# To avoid this, it is crucial to use a validation dataset for model selection and hyperparameter tuning, 
# leaving the test set untouched until the very end as a final evaluation of the model.



