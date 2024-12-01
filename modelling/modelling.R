library(dplyr)
library(reshape2)
library(tidyr)
library(ranger)
library(caret)
library(keras)
library(tensorflow)
library(iml)
library(ggplot2)

setwd("/Users/oliver/RStudio")

# Function to perform one iteration of Random Forest with training and evaluation
perform_iteration <- function(data, response_var, predictor_vars, num_trees, folds) {
  # Sample the data
  sampled_data <- data %>% sample_frac(sample_fraction)
  
  # Split the data into training and testing sets
  train_index <- createDataPartition(sampled_data[[response_var]], p = 0.8, list = FALSE)
  train_data <- sampled_data[train_index, ]
  test_data <- sampled_data[-train_index, ]
  
  # Set up cross-validation controls
  train_control <- trainControl(method = "cv", number = folds)
  
  # Create tuning grid for hyperparameters
  tune_grid <- expand.grid(
    mtry = floor(sqrt(length(predictor_vars))),
    splitrule = "variance",
    min.node.size = 1
  )
  
  # Train the Random Forest model
  ranger_model <- caret::train(
    as.formula(paste(response_var, "~ .")),
    data = train_data,
    method = "ranger",
    trControl = train_control,
    tuneGrid = tune_grid,
    num.trees = num_trees,
    importance = 'impurity'
  )
  
  # Make predictions on the test data
  predictions <- predict(ranger_model, newdata = test_data)
  test_rmse <- sqrt(mean((predictions - test_data[[response_var]])^2))
  
  # Calculate R-squared
  ss_total <- sum((test_data[[response_var]] - mean(test_data[[response_var]]))^2)
  ss_residual <- sum((predictions - test_data[[response_var]])^2)
  r_squared <- 1 - (ss_residual / ss_total)
  
  list(model = ranger_model,
       test_rmse = test_rmse,
       r_squared = r_squared,
       importance = ranger_model$finalModel$variable.importance)
}

# Load and prepare the data
df <- read.delim("final_quantum.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Set sgRNAID as row names and remove it from the data frame
rownames(df) <- df$sgRNAID
df$sgRNAID <- NULL  

# Define the response variable and predictor variables
response_var <- 'cut.score'
predictor_vars <- setdiff(names(df), response_var)

# Check the Pearson correlation
# Calculate the correlation matrix
correlation_matrix <- cor(df[, predictor_vars], use = "complete.obs")

# Find pairs with Pearson correlation larger than 0.9
high_correlation_pairs <- which(abs(correlation_matrix) > 0.9 & abs(correlation_matrix) < 1, arr.ind = TRUE)

# Remove duplicate pairs and extract unique variable names
high_correlation_pairs <- unique(t(apply(high_correlation_pairs, 1, sort)))
high_correlation_vars <- unique(c(rownames(high_correlation_pairs), colnames(high_correlation_pairs)))

# Print variables and their correlation values
for (i in seq_len(nrow(high_correlation_pairs))) {
  var1 <- rownames(correlation_matrix)[high_correlation_pairs[i, 1]]
  var2 <- colnames(correlation_matrix)[high_correlation_pairs[i, 2]]
  corr_value <- correlation_matrix[var1, var2]
  
  cat(sprintf("Variables: %s and %s - Pearson correlation: %.2f\n", var1, var2, corr_value))
}

# Parameters for the iterative random forest
iterations <- 10
num_trees <- 400
folds <- 5
sample_fraction <- 0.2

# Initialize lists to store results
importance_list <- vector("list", iterations)
test_rmse_list <- numeric(iterations)
r_squared_list <- numeric(iterations)
adjusted_r_squared_list <- numeric(iterations)

# Initialize variables to store the model after the final iteration
final_ranger_model <- NULL

# Convert predictor variables to numeric
df[predictor_vars] <- lapply(df[predictor_vars], as.numeric)

# Perform iterations
for (i in 1:iterations) {
  iteration_result <- perform_iteration(df, response_var, predictor_vars, num_trees, folds)
  importance_list[[i]] <- iteration_result$importance
  test_rmse_list[i] <- iteration_result$test_rmse
  r_squared_list[i] <- iteration_result$r_squared
  
  # Save the model from the current iteration
  final_ranger_model <- iteration_result$model
  
  # Calculate the adjusted R-squared
  n <- nrow(df)
  p <- length(predictor_vars)
  adjusted_r_squared <- 1 - ((1 - iteration_result$r_squared) * (n - 1) / (n - p - 1))
  adjusted_r_squared_list[i] <- adjusted_r_squared
}

# Save the iterative Random Forest model
saveRDS(final_ranger_model, file = "iRF.rds")

# Print mean test RMSE, R-squared
mean_test_rmse <- mean(test_rmse_list)
mean_r_squared <- mean(r_squared_list)
mean_adjusted_r_squared <- mean(adjusted_r_squared_list)
cat("Mean Test RMSE:", mean_test_rmse, "\n")
cat("Mean R-squared:", mean_r_squared, "\n")
cat("Mean Adjusted R-squared:", mean_adjusted_r_squared, "\n")

# Plot RMSE across iterations
plot(test_rmse_list, type = "b", main = "RMSE Across Iterations", xlab = "Iteration", ylab = "RMSE")

# Plot adjusted R-squared across iterations
plot(adjusted_r_squared_list, type = "b", main = "Adjusted R-squared Across Iterations", xlab = "Iteration", ylab = "adjusted R-squared")







##### Evaluate model and feature importance
# Calculate mean feature importance across iterations
mean_importance <- Reduce("+", importance_list) / length(importance_list)

# Save the mean importance to a CSV file
importance_df <- data.frame(feature = names(mean_importance), importance = mean_importance)
write.csv(importance_df, "feature_importance.csv", row.names = FALSE)

# Combine the importance with signs
importance_with_sign <- data.frame(
  feature = names(mean_importance),
  importance = mean_importance,
  sign = sign(cor(df[,predictor_vars], df$cut.score))
)

# Sort by absolute importance and select the top 50
top_50_features <- importance_with_sign %>%
  mutate(importance = importance * sign) %>%
  arrange(desc(abs(importance))) %>%
  head(50)

# Sort by absolute importance and select the least important 50
least_50_features <- importance_with_sign %>%
  mutate(importance = importance * sign) %>%
  arrange(abs(importance)) %>%
  head(50)

# Determine the direction of effect using correlation
direction_of_effect <- sapply(predictor_vars, function(var) {
  cor(df[[var]], df[[response_var]], use = "complete.obs")
})

# Add direction to the importance dataframe
importance_df <- importance_df %>%
  mutate(direction = ifelse(direction_of_effect[feature] > 0, "pos", "neg"))

# Normalize the importance
top_50_features <- top_50_features %>%
  mutate(normalized_importance = importance / max(abs(importance)))

least_50_features <- least_50_features %>%
  mutate(normalized_importance = importance / max(abs(importance)))

# Plot top 50 features
ggplot(top_50_features, aes(x = reorder(feature, normalized_importance), y = normalized_importance, fill = factor(sign))) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(
    title = "Top 50 Feature Importance",
    x = "Feature",
    y = "Normalized Importance",
    fill = "Direction of Effect"
  ) +
  scale_fill_manual(values = c("-1" = "red", "1" = "blue"), labels = c("Negative", "Positive")) +
  theme_minimal()
      

# Plot least important 50 features
ggplot(least_50_features, aes(x = reorder(feature, normalized_importance), y = normalized_importance, fill = factor(sign))) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(
    title = "Least 50 Feature Importance",
    x = "Feature",
    y = "Normalized Importance",
    fill = "Direction of Effect"
  ) +
  scale_fill_manual(values = c("-1" = "red", "1" = "blue"), labels = c("Negative", "Positive")) +
  theme_minimal()








###### Evaluate the effect of adding quantum chemical features
# Load and prepare the raw.matrix.txt data
df_raw <- read.delim("raw.matrix.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Set sgRNAID as row names
rownames(df_raw) <- df_raw$sgRNAID
df_raw$sgRNAID <- NULL  # Remove sgRNAID column

# Define the response variable and predictor variables
response_var_raw <- 'cut.score'
predictor_vars_raw <- setdiff(names(df_raw), response_var_raw)  # Exclude response_var from predictors

# Convert predictor variables to numeric
df_raw[predictor_vars_raw] <- lapply(df_raw[predictor_vars_raw], as.numeric)

# Parameters for the iterative random forest
iterations <- 10
num_trees <- 400
folds <- 5
sample_fraction <- 0.2

# Initialize lists to store results
importance_list_raw <- vector("list", iterations)
test_rmse_list_raw <- numeric(iterations)
r_squared_list_raw <- numeric(iterations)
adjusted_r_squared_list_raw <- numeric(iterations)

perform_iteration <- function(data, response_var, predictor_vars, num_trees, folds) {
  sampled_data <- data %>% sample_frac(sample_fraction)
  
  train_index <- createDataPartition(sampled_data[[response_var]], p = 0.8, list = FALSE)
  train_data <- sampled_data[train_index, ]
  test_data <- sampled_data[-train_index, ]
  
  train_control <- trainControl(method = "cv", number = folds)
  
  tune_grid <- expand.grid(
    mtry = floor(sqrt(length(predictor_vars))),
    splitrule = "variance",
    min.node.size = 1
  )
  
  ranger_model <- caret::train(
    as.formula(paste(response_var, "~ .")),
    data = train_data,
    method = "ranger",
    trControl = train_control,
    tuneGrid = tune_grid,
    num.trees = num_trees,
    importance = 'impurity'
  )
  
  predictions <- predict(ranger_model, newdata = test_data)
  test_rmse <- sqrt(mean((predictions - test_data[[response_var]])^2))
  
  # Calculate R-squared
  ss_total <- sum((test_data[[response_var]] - mean(test_data[[response_var]]))^2)
  ss_residual <- sum((predictions - test_data[[response_var]])^2)
  r_squared <- 1 - (ss_residual / ss_total)
  
  list(model = ranger_model,
       test_rmse = test_rmse,
       r_squared = r_squared,
       importance = ranger_model$finalModel$variable.importance)
}

# Perform iterations
for (i in 1:iterations) {
  iteration_result <- perform_iteration(df_raw, response_var_raw, predictor_vars_raw, num_trees, folds)
  importance_list_raw[[i]] <- iteration_result$importance
  test_rmse_list_raw[i] <- iteration_result$test_rmse
  r_squared_list_raw[i] <- iteration_result$r_squared
  
  # Calculate the adjusted R-squared
  n <- nrow(df_raw)
  p <- length(predictor_vars)
  adjusted_r_squared <- 1 - ((1 - iteration_result$r_squared) * (n - 1) / (n - p - 1))
  adjusted_r_squared_list[i] <- adjusted_r_squared
}

# Print mean test RMSE and R-squared for conventional model
mean_test_rmse_raw <- mean(test_rmse_list_raw)
mean_r_squared_raw <- mean(r_squared_list_raw)
mean_adjusted_r_squared_raw <- mean(adjusted_r_squared_list)
cat("Mean Test RMSE (conventional):", mean_test_rmse_raw, "\n")
cat("Mean R-squared (conventional):", mean_r_squared_raw, "\n")
cat("Mean Adjusted R-squared (conventional):", mean_adjusted_r_squared_raw, "\n")

cat("Comparison of RMSE:")
cat("Conventional Model RMSE:", mean_test_rmse_raw, "\n")
cat("Quantum Model RMSE:", mean_test_rmse, "\n")

cat("Comparison of Adjusted R-squared:")
cat("Conventional Model Adjusted R-squared:", mean_adjusted_r_squared_raw, "\n")
cat("Quantum Model Adjusted R-squared:", mean_adjusted_r_squared, "\n")


# Calculate mean feature importance for the conventional model
mean_importance_raw <- Reduce("+", importance_list_raw) / length(importance_list_raw)
importance_df_raw <- data.frame(feature = names(mean_importance_raw), importance = mean_importance_raw)

# Determine the direction of effect
direction_of_effect <- sapply(predictor_vars, function(var) {
  cor(df[[var]], df[[response_var]], use = "complete.obs")
})

direction_of_effect <- direction_of_effect[importance_df_raw$feature]

# Combine the importance with signs and directions
importance_with_sign_raw <- importance_df_raw %>%
  mutate(sign = sign(direction_of_effect)) %>%
  mutate(direction = ifelse(direction_of_effect > 0, "pos", "neg"))

# Sort by absolute importance and select the top 50
top_50_features_raw <- importance_with_sign_raw %>%
  mutate(importance = importance * sign) %>%
  arrange(desc(abs(importance))) %>%
  head(50)

# Normalize the importance
top_50_features_raw <- top_50_features_raw %>%
  mutate(normalized_importance = importance / max(abs(importance)))

# Plot the top 50 features
ggplot(top_50_features_raw, aes(x = reorder(feature, normalized_importance), y = normalized_importance, fill = factor(sign))) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(
    title = "Top 50 Feature Importance (Conventional Model)",
    x = "Feature",
    y = "Normalized Importance",
    fill = "Direction of Effect"
  ) +
  scale_fill_manual(values = c("-1" = "red", "1" = "blue"), labels = c("Negative", "Positive")) +
  theme_minimal()










##### Training of a neural network model
# Select top 100 features
top_100_features <- importance_df %>%
  arrange(desc(importance)) %>%
  head(100) %>%
  pull(feature)

df_top100 <- df[, c("cut.score", top_100_features)]

# Normalize the data
preproc <- preProcess(df_top100[, -1], method = c("center", "scale"))
df_top100_normalized <- predict(preproc, df_top100)

# Separate the target variable and features
x <- as.matrix(df_top100_normalized[, -1])
y <- as.numeric(df_top100_normalized$cut.score)

# Cross-validation parameters
k <- 5  # Number of folds
folds <- cut(seq(1, nrow(df_top100_normalized)), breaks = k, labels = FALSE)

# Initialize vectors to store metrics
cv_mse <- c()
cv_mae <- c()
cv_r2 <- c()
cv_adjusted_r2 <- c()

# Perform cross-validation
for (i in 1:k) {
  # Split the data into training and validation sets
  validation_indices <- which(folds == i, arr.ind = TRUE)
  x_train <- x[-validation_indices,]
  y_train <- y[-validation_indices]
  x_val <- x[validation_indices,]
  y_val <- y[validation_indices]
  
  # Define the neural network model
  model <- keras_model_sequential() %>%
    # First layer
    layer_dense(units = 64, input_shape = ncol(x), 
                kernel_regularizer = regularizer_l2(0.001)) %>%
    layer_batch_normalization() %>%
    layer_activation('relu') %>%
    layer_dropout(rate = 0.5) %>%
    
    # Second layer
    layer_dense(units = 32, 
                kernel_regularizer = regularizer_l2(0.001)) %>%
    layer_batch_normalization() %>%
    layer_activation('relu') %>%
    layer_dropout(rate = 0.5) %>%
    
    # Third layer
    layer_dense(units = 16, 
                kernel_regularizer = regularizer_l2(0.001)) %>%
    layer_batch_normalization() %>%
    layer_activation('relu') %>%
    layer_dropout(rate = 0.2) %>%
    
    # Output layer
    layer_dense(units = 1)
  
  # Compile the model
  model %>% compile(
    loss = 'mean_squared_error',
    optimizer = tf$keras$optimizers$legacy$Adam(learning_rate = 0.001),
    metrics = c('mean_absolute_error')
  )
  
  # Train the model with early stopping
  callback_early_stopping <- callback_early_stopping(monitor = 'val_loss', patience = 10)
  
  history <- model %>% fit(
    x_train, y_train,
    epochs = 100,
    batch_size = 32,
    validation_data = list(x_val, y_val),
    callbacks = list(callback_early_stopping)
  )
  
  # Calculate valuation metrics
  val_predictions <- model %>% predict(x_val)
  mse <- mean((val_predictions - y_val)^2)
  mae <- mean(abs(val_predictions - y_val))
  sst <- sum((y_val - mean(y_val))^2)
  ssr <- sum((y_val - val_predictions)^2)
  r_squared <- 1 - (ssr / sst)

  n <- length(y_val)
  p <- ncol(x_train)
  adjusted_r_squared <- 1 - ((1 - r_squared) * (n - 1) / (n - p - 1))
  
  # Store the metrics
  cv_mse <- c(cv_mse, mse)
  cv_mae <- c(cv_mae, mae)
  cv_r2 <- c(cv_r2, r_squared)
  cv_adjusted_r2 <- c(cv_adjusted_r2, adjusted_r_squared)
}

# Report cross-validated metrics
cat("Cross-validated MSE:", mean(cv_mse), "\n")
cat("Cross-validated MAE:", mean(cv_mae), "\n")
cat("Cross-validated R-squared:", mean(cv_r2), "\n")
cat("Cross-validated Adjusted R-squared:", mean(cv_adjusted_r2), "\n")


# Final model training on the full dataset
final_model <- keras_model_sequential() %>%
  # First layer
  layer_dense(units = 64, input_shape = ncol(x), 
              kernel_regularizer = regularizer_l2(0.001)) %>%
  layer_batch_normalization() %>%
  layer_activation('relu') %>%
  layer_dropout(rate = 0.5) %>%
  
  # Second layer
  layer_dense(units = 32, 
              kernel_regularizer = regularizer_l2(0.001)) %>%
  layer_batch_normalization() %>%
  layer_activation('relu') %>%
  layer_dropout(rate = 0.5) %>%
  
  # Third layer
  layer_dense(units = 16, 
              kernel_regularizer = regularizer_l2(0.001)) %>%
  layer_batch_normalization() %>%
  layer_activation('relu') %>%
  layer_dropout(rate = 0.2) %>%
  
  # Output layer
  layer_dense(units = 1)

final_model %>% compile(
  loss = 'mean_squared_error',
  optimizer = tf$keras$optimizers$legacy$Adam(learning_rate = 0.001),
  metrics = c('mean_absolute_error')
)

final_history <- final_model %>% fit(
  x, y,
  epochs = 100,
  batch_size = 64,
  validation_split = 0.2,
  callbacks = list(callback_early_stopping)
)

# Plot training and validation metrics
plot(final_history)

# Calculate valuation metrics
final_predictions <- final_model %>% predict(x)
final_mse <- mean((final_predictions - y)^2)
final_rmse <- sqrt(final_mse)
final_mae <- mean(abs(final_predictions - y))
final_sst <- sum((y - mean(y))^2)
final_ssr <- sum((y - final_predictions)^2)
final_r_squared <- 1 - (final_ssr / final_sst)

n_final <- length(y)
p_final <- ncol(x)
final_adjusted_r_squared <- 1 - ((1 - final_r_squared) * (n_final - 1) / (n_final - p_final - 1))

# Print final evaluation metrics
cat("Final MSE:", final_mse, "\n")
cat("Final RMSE:", final_rmse, "\n")
cat("Final MAE:", final_mae, "\n")
cat("Final R-squared:", final_r_squared, "\n")
cat("Final Adjusted R-squared:", final_adjusted_r_squared, "\n")



###### Calculate feature importance with SHAP value
# Convert trained model to a predictor object
model_predictor <- Predictor$new(final_model, data = as.data.frame(x), y = y)

# Calculate SHAP values
correct_x <- as.data.frame(x[1, , drop = FALSE])
shap_values <- Shapley$new(model_predictor, x.interest = correct_x)

shap_df <- as.data.frame(shap_values$results)
shap_df$phi <- as.numeric(shap_df$phi)

# Sort the SHAP values
shap_df <- shap_df[order(abs(shap_df$phi), decreasing = TRUE), ]

# Select the top 50 SHAP values
top_50_shap <- shap_df[1:50, ]

# Plot top 50 SHAP values
ggplot(top_50_shap, aes(x = reorder(feature, phi), y = phi, fill = phi > 0)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(
    title = "Top 50 SHAP Values",
    x = "Feature",
    y = "SHAP Value"
  ) +
  scale_fill_manual(values = c("TRUE" = "blue", "FALSE" = "red"), labels = c("Negative", "Positive")) +
  theme_minimal()

# Save the Neural Network model
save_model_tf(final_model, filepath = "neural_network.keras")








###### Validation on human dataset
# Load the Random Forest model
loaded_ranger_model <- readRDS("iRF.rds")

# Load the Neural Network model
loaded_nn_model <- load_model_tf("neural_network.keras")

# Load the validation dataset
val_data <- read.delim("val_final_quantum.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Set sgRNAID as row names and remove it from the data frame
rownames(val_data) <- val_data$sgRNAID
val_data$sgRNAID <- NULL  # Remove sgRNAID column

# Ensure predictor variables are numeric
predictor_vars <- setdiff(names(val_data), "cut.score")
val_data[predictor_vars] <- lapply(val_data[predictor_vars], as.numeric)

# Separate features and target variable
x_val <- val_data[, predictor_vars]
x_val_selected <- x_val[, top_100_features]
y_val <- val_data$cut.score

# Number of observations and predictors
n <- nrow(val_data)
p <- length(predictor_vars)

# Normalize the validation data
x_val_normalized <- predict(preproc, x_val)
x_val_selected_normalized <- predict(preproc, x_val_selected)

# Make predictions with the Random Forest model
rf_predictions <- predict(loaded_ranger_model, newdata = x_val_normalized)

# Make predictions with the Neural Network model
nn_predictions <- loaded_nn_model %>% predict(as.matrix(x_val_selected_normalized))

# Calculate evaluation metrics for Random Forest
rf_mse <- mean((rf_predictions - y_val)^2)
rf_rmse <- sqrt(rf_mse)
rf_mae <- mean(abs(rf_predictions - y_val))
rf_r_squared <- 1 - (sum((y_val - rf_predictions)^2) / sum((y_val - mean(y_val))^2))
rf_adjusted_r_squared <- 1 - ((1 - rf_r_squared) * (n - 1) / (n - p - 1))

cat("Random Forest Model Evaluation:\n")
cat("RMSE:", rf_rmse, "\n")
cat("MAE:", rf_mae, "\n")
cat("R-squared:", rf_r_squared, "\n")
cat("Adjusted R-squared:", rf_adjusted_r_squared, "\n")

# Calculate evaluation metrics for Neural Network
nn_mse <- mean((nn_predictions - y_val)^2)
nn_rmse <- sqrt(nn_mse)
nn_mae <- mean(abs(nn_predictions - y_val))
nn_r_squared <- 1 - (sum((y_val - nn_predictions)^2) / sum((y_val - mean(y_val))^2))
nn_adjusted_r_squared <- 1 - ((1 - nn_r_squared) * (n - 1) / (n - p - 1))

cat("Neural Network Model Evaluation:\n")
cat("RMSE:", nn_rmse, "\n")
cat("MAE:", nn_mae, "\n")
cat("R-squared:", nn_r_squared, "\n")
cat("Adjusted R-squared:", nn_adjusted_r_squared, "\n")

