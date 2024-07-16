library(BART)
library(survival)
library(survminer)
library(ggplot2)
library(rpart)
library(rpart.plot)







##### First model (complete remission) #####

# Select columns excluding 'complete_remission', 'as.factor.myRsp.', 'X', and 'bcr_patient_barcode'
pred1 <- data1[, !colnames(data1) %in% c("complete_remission", "treatment_outcome_first_course", "as.factor.myRsp.", "X", "bcr_patient_barcode", "partial_remission")]

# Ensure pred1 is a matrix
pred1 <- as.matrix(pred1)

# # Fit the BART
bart_model_complete <- pbart(x.train = pred1, y.train = data1$complete_remission,
                              ndpost = 5000, ntree = 1, base=.5)
var_complete <- bart_model_complete$varcount.mean

# Fit the CART model
cart_model_complete <- rpart(complete_remission ~ ., data = data.frame(complete_remission = data1$complete_remission, pred1), method = "class")
summary(cart_model_complete)
rpart.plot(cart_model_complete, main = "Complete remission vs. rest (without interaction terms)")



##### Second model (partial remission) #####

# Select columns excluding 'complete_remission', 'as.factor.myRsp.', 'X', and 'bcr_patient_barcode'
pred2 <- data2[, !colnames(data2) %in% c("complete_remission", "treatment_outcome_first_course", "as.factor.myRsp.", "X", "bcr_patient_barcode", "partial_remission")]

# Ensure pred1 is a matrix
pred2 <- as.matrix(pred2)

# Fit the BART
bart_model_partial <- pbart(x.train = pred2, y.train = data2$partial_remission,
                             ndpost = 5000, ntree = 1, base = .80)
var_partial <- bart_model_partial$varcount.mean
barplot(var_partial)

treeDepth(bart_model_partial)

# Fit the CART model
cart_model_partial <- rpart(partial_remission ~ ., data = data.frame(partial_remission = data2$partial_remission, pred2), method = "class")
summary(cart_model_partial)
rpart.plot(cart_model_partial, main = "Partial remission vs. rest (without interaction terms)")



##### Third model (stable disease) #####

# Select columns excluding 'complete_remission', 'as.factor.myRsp.', 'X', and 'bcr_patient_barcode'
pred3 <- data3[, !colnames(data3) %in% c("complete_remission", "treatment_outcome_first_course", "as.factor.myRsp.", "X", "bcr_patient_barcode","stable_disease", "partial_remission")]

# Ensure pred1 is a matrix
pred3 <- as.matrix(pred3)

# Fit the BART
bart_model_stable <- pbart(x.train = pred3, y.train = data3$stable_disease,
                            ndpost = 5000, ntree = 1, base = 0.95)
var_stable <- bart_model_stable$varcount.mean
barplot(var_stable)

cart_model_stable <- rpart(stable_disease ~ ., data = data.frame(stable_disease = data3$stable_disease, pred3), method = "class")
summary(cart_model_stable)
rpart.plot(cart_model_stable, main = "Stable disease vs. rest (without interaction terms)")




# make df with new interaction colums
# Function to create interaction terms and combine them with the original matrix
create_interaction_matrix <- function(data, interaction_col) {
  # Convert the data to a data frame
  data_df <- as.data.frame(data)
  
  # Identify all protein columns (excluding the interaction column and other non-protein columns)
  protein_cols <- setdiff(colnames(data_df), interaction_col)
  
  # Create interaction terms between the interaction column and each protein column
  interaction_terms <- sapply(protein_cols, function(protein) {
    interaction_term <- data_df[[interaction_col]] * data_df[[protein]]
    colname <- paste(interaction_col, protein, sep = "_x_")
    interaction_term
  })
  
  # Convert interaction terms to a data frame and name the columns
  interaction_df <- as.data.frame(interaction_terms)
  colnames(interaction_df) <- paste(interaction_col, protein_cols, sep = "_x_")
  
  # Combine the original matrix with the new interaction terms
  new_data <- cbind(data_df, interaction_df)
  
  # Convert to matrix
  new_data_matrix <- as.matrix(new_data)
  
  return(new_data_matrix)
}

# Apply the function to pred1 and pred3
pred1 <- as.matrix(pred1)
pred2 <- as.matrix(pred2)
pred3 <- as.matrix(pred3)

new_pred1 <- create_interaction_matrix(pred1, "newTRT")
new_pred2 <- create_interaction_matrix(pred2, "newTRT")
new_pred3 <- create_interaction_matrix(pred3, "newTRT")


# Fit the new second CART model
cart_model_partial <- rpart(partial_remission ~ ., data = data.frame(partial_remission = data2$partial_remission, new_pred2), method = "class")
summary(cart_model_partial)
rpart.plot(cart_model_partial, main = "Partial remission vs. rest (w/ interaction terms)")

# Fit the new first CART model
cart_model_complete <- rpart(complete_remission ~ ., data = data.frame(complete_remission = data1$complete_remission, new_pred1), method = "class")
summary(cart_model_complete)
rpart.plot(cart_model_complete, main = "Complete remission vs. rest (w/ interaction terms)")

# Fit the new third CART model
cart_model_stable <- rpart(stable_disease ~ ., data = data.frame(stable_disease = data3$stable_disease, new_pred3), method = "class")
summary(cart_model_stable)
rpart.plot(cart_model_stable, main = "Stable vs. rest (w/ interaction terms)")













# Install and load necessary packages
if (!require("bartMachine")) install.packages("bartMachine")
library(bartMachine)

# Assuming pred1 and data1$complete_remission are already defined and prepared

# Fit the BART model
set_bart_machine_num_cores(1) # Set number of cores to 1 for reproducibility
bart_model <- bartMachine(X = pred1, y = data1$complete_remission, num_trees = 1, num_burn_in = 1000, num_iterations_after_burn_in = 5000)

# Print the summary of the model
print(bart_model)

# Plot the variable importance
investigate_var_importance(bart_model)

# If you want to plot the tree, you will need to use internal functions
# Note: BART models are ensemble models, and plotting a single tree might not represent the entire model well.

# Extract and plot a single tree (for the first iteration)
single_tree <- bart_model$java_bart_machine$bart_trees$`0`

# Plot the single tree using a simple plot function
plot(single_tree)
text(single_tree, use.n = TRUE, all = TRUE, cex = 0.8)











# Install and load necessary packages
if (!require("rpart")) install.packages("rpart")
library(rpart)

if (!require("rpart.plot")) install.packages("rpart.plot")
library(rpart.plot)

# Assuming pred1 is your existing matrix and data1 is the data frame containing 'newTRT'
# Prepare the data by creating a data frame
data_combined <- data.frame(newTRT = data1$newTRT, pred1)

# Fit the CART model
cart_model_confounders <- rpart(newTRT ~ ., data = data.frame(pred1), method = "class")

# Print the summary of the model
summary(cart_model_confounders)

# Plot the tree
plot(cart_model_confounders)
text(cart_model_confounders, use.n = TRUE, all = TRUE, cex = 0.8)

# Optional: Better plot with rpart.plot
rpart.plot(cart_model_confounders, type = 2, extra = 104, under = TRUE, faclen = 0)

