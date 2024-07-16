library(dplyr)
library(BART)
library(rpart)
library(rpart.plot)
library(corrplot)
library(Matrix)

data1 <- read.csv("~/Desktop/Amsterdam UMC/Data Research Camp 2024/data1.csv")[,-(1:2)]
data2 <- read.csv("~/Desktop/Amsterdam UMC/Data Research Camp 2024/data2.csv")[,-(1:2)]
data3 <- read.csv("~/Desktop/Amsterdam UMC/Data Research Camp 2024/data3.csv")[,-(1:2)]


# Function to create interaction matrix
create_interaction_matrix <- function(data, interaction_col) {
  data_df <- as.data.frame(data)
  protein_cols <- setdiff(colnames(data_df), interaction_col)
  interaction_terms <- sapply(protein_cols, function(protein) {
    interaction_term <- data_df[[interaction_col]] * data_df[[protein]]
    colname <- paste(interaction_col, protein, sep = "_x_")
    interaction_term
  })
  interaction_df <- as.data.frame(interaction_terms)
  colnames(interaction_df) <- paste(interaction_col, protein_cols, sep = "_x_")
  new_data <- cbind(data_df, interaction_df)
  new_data_matrix <- as.matrix(new_data)
  return(new_data_matrix)
}

# Loading data
data1 <- read.csv("~/Desktop/Amsterdam UMC/Data Research Camp 2024/data1.csv")[,-(1:2)]
data2 <- read.csv("~/Desktop/Amsterdam UMC/Data Research Camp 2024/data2.csv")[,-(1:2)]
data3 <- read.csv("~/Desktop/Amsterdam UMC/Data Research Camp 2024/data3.csv")[,-(1:2)]

# Columns
treat_col <- c("radiation_treatment_adjuvant", "targeted_molecular_therapy", "both_treatments")
protein_col <- c("ACVRL1.R.C", "Src_pY416.R.C", "HSP70.R.C", "HER2_pY1248.R.C", "PAI.1.M.E",
                 "Smad1.R.V", "FOXO3a_pS318_S321.R.C", "PRAS40_pT246.R.V", "Cyclin_E2.R.C",
                 "SF2.M.V", "Bad_pS112.R.V", "Lck.R.V", "Caspase.7_cleavedD198.R.C",
                 "Paxillin.R.C", "MYH11.R.V", "X14.3.3_epsilon.M.C", "Akt_pS473.R.V",
                 "Caveolin.1.R.V", "Rab25.R.V", "YAP_pS127.R.E", "RBM15.R.V", "Claudin.7.R.V",
                 "ER.alpha.R.V", "C.Raf.R.V", "CD31.M.V", "Ku80.R.C", "Bcl.2.M.V", "GSK3.alpha.beta.M.V")

# Function to fit BART and CART models
fit_models <- function(data, response_col, predictors, interaction_col = FALSE) {
  data_tree <- data[, predictors]
  response <- data[[response_col]]
  
  # Fit BART model
  set.seed(1)
  bart_model <- pbart(x.train = data_tree, y.train = response, ndpost = 5000, ntree = 1, base = 0.9)
  var_stable <- bart_model$varcount.mean
  barplot(var_stable, main = paste(response_col, "vs. else"))
  
  # Fit CART model
  data_cart <- data.frame(data_tree, response = response)
  cart_model <- rpart(response ~ ., data = data_cart, method = "class")
  rpart.plot(cart_model, main = paste("CART for", response_col))
  
  # Variable importance from CART model
  var_importance <- cart_model$variable.importance
  barplot(var_importance, main = paste("Variable Importance (CART) for", response_col), col = "blue", las = 2)
  
  # Fit BART model with interactions if interaction_col is provided
  if (interaction_col) {
    data_tree_int <- create_interaction_matrix(data_tree, interaction_col)
    bart_model_int <- pbart(x.train = data_tree_int, y.train = response, ndpost = 5000, ntree = 1, base = 0.9)
    var_stable_int <- bart_model_int$varcount.mean
    barplot(var_stable_int, main = paste(response_col, "vs. else (with interactions)"))
  }
}

# Progressive vs else (WITHOUT interaction terms)
fit_models(data1, 
           "progressive_disease", 
           colnames(data1)[!colnames(data1) %in% c("treatment_outcome_first_course", "as.factor.myRsp.", "radiation_treatment_adjuvant", "targeted_molecular_therapy", "both_treatments", "progressive_disease")],
           interaction_col = FALSE)

# Progressive vs else (WITH interaction terms)
fit_models(data1, 
           "progressive_disease", 
           colnames(data1)[!colnames(data1) %in% c("treatment_outcome_first_course", "as.factor.myRsp.", "radiation_treatment_adjuvant", "targeted_molecular_therapy", "both_treatments", "progressive_disease")],
           interaction_col = TRUE)

# Stable disease vs else (WITHOUT interaction terms)
fit_models(data2, 
           "stable_disease", 
           colnames(data2)[!colnames(data2) %in% c("treatment_outcome_first_course", "as.factor.myRsp.", "radiation_treatment_adjuvant", "targeted_molecular_therapy", "both_treatments", "progressive_disease", "stable_disease")],
           interaction_col = FALSE)

# Stable disease vs else (WITH interaction terms)
fit_models(data2, 
           "stable_disease", 
           colnames(data2)[!colnames(data2) %in% c("treatment_outcome_first_course", "as.factor.myRsp.", "radiation_treatment_adjuvant", "targeted_molecular_therapy", "both_treatments", "progressive_disease", "stable_disease")],
           interaction_col = TRUE)

# Partial remission vs else (WITHOUT interaction terms)
fit_models(data3, "partial_remission", colnames(data3)[!colnames(data3) %in% c("treatment_outcome_first_course", "as.factor.myRsp.", "radiation_treatment_adjuvant", "targeted_molecular_therapy", "both_treatments", "progressive_disease", "stable_disease", "partial_remission")],
           interaction_col = FALSE)

# Partial remission vs else (WITH interaction terms)
fit_models(data3, "partial_remission", colnames(data3)[!colnames(data3) %in% c("treatment_outcome_first_course", "as.factor.myRsp.", "radiation_treatment_adjuvant", "targeted_molecular_therapy", "both_treatments", "progressive_disease", "stable_disease", "partial_remission")],
           interaction_col = TRUE)
