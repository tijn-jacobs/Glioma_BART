## Research Data Camp
setwd("~/Desktop/RDC")
rm(list = ls())
library(dplyr)
library(BART)
library(rpart)
library(rpart.plot)
library(corrplot)
library(Matrix)

data1 <- read.csv("data1.csv")[,-(1:2)]
data2 <- read.csv("data2.csv")[,-(1:2)]
data3 <- read.csv("data3.csv")[,-(1:2)]

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

# Function to fit BART and CART models
fit_models <- function(data, ndpost, ntree, seed, response_col, predictors, interaction_col = FALSE) {
  data_tree <- data[, predictors]
  response <- data[[response_col]]
  
  # Fit BART model
  set.seed(seed)
  bart_model <- pbart(x.train = data_tree, y.train = response, ndpost = ndpost, ntree = ntree, base = 0.9)
  var_stable <- bart_model$varcount.mean
  barplot(var_stable, main = paste(response_col, "vs. else"))
  
  # Fit CART model
  data_cart <- data.frame(data_tree, response = response)
  cart_model <- rpart(response ~ ., data = data_cart, method = "class")
  rpart.plot(cart_model, main = paste("CART for", response_col))
  
  # Variable importance from CART model
  var_importance <- cart_model$variable.importance
  barplot(var_importance, main = paste("Variable Importance (CART) for", response_col), col = "blue", las = 2)
  
  # Fit BART and CART models with interactions if interaction_col is provided
  if (interaction_col != FALSE) {
    data_tree_int <- create_interaction_matrix(data_tree, interaction_col)
    bart_model_int <- pbart(x.train = data_tree_int, y.train = response, ndpost = ndpost, ntree = ntree, base = 0.9)
    var_stable_int <- bart_model_int$varcount.mean
    barplot(var_stable_int, main = paste(response_col, "vs. else (with interactions)"))
    
    # Fit CART model with interactions
    data_cart_int <- data.frame(data_tree_int, response = response)
    cart_model_int <- rpart(response ~ ., data = data_cart_int, method = "class")
    rpart.plot(cart_model_int, main = paste("CART with Interactions for", response_col))
    
    # Variable importance from CART model with interactions
    var_importance_int <- cart_model_int$variable.importance
    barplot(var_importance_int, main = paste("Variable Importance (CART with Interactions) for", response_col), col = "blue", las = 2)
  }
}

##
# Columns
treat_col <- c("radiation_treatment_adjuvant", "targeted_molecular_therapy", "both_treatments")
protein_col <- c("ACVRL1.R.C", "Src_pY416.R.C", "HSP70.R.C", "HER2_pY1248.R.C", "PAI.1.M.E",
                 "Smad1.R.V", "FOXO3a_pS318_S321.R.C", "PRAS40_pT246.R.V", "Cyclin_E2.R.C",
                 "SF2.M.V", "Bad_pS112.R.V", "Lck.R.V", "Caspase.7_cleavedD198.R.C",
                 "Paxillin.R.C", "MYH11.R.V", "X14.3.3_epsilon.M.C", "Akt_pS473.R.V",
                 "Caveolin.1.R.V", "Rab25.R.V", "YAP_pS127.R.E", "RBM15.R.V", "Claudin.7.R.V",
                 "ER.alpha.R.V", "C.Raf.R.V", "CD31.M.V", "Ku80.R.C", "Bcl.2.M.V", "GSK3.alpha.beta.M.V")

# Progressive vs else (WITHOUT interaction terms)
fit_models(data1, ndpost = 5000, ntree = 1, seed = 42, 
           "progressive_disease", 
           colnames(data1)[!colnames(data1) %in% c("treatment_outcome_first_course", "as.factor.myRsp.", "radiation_treatment_adjuvant", "targeted_molecular_therapy", "both_treatments", "progressive_disease")],
           interaction_col = FALSE)

# Progressive vs else (WITH interaction terms)
fit_models(data1, ndpost = 5000, ntree = 1, seed = 42, 
           "progressive_disease", 
           colnames(data1)[!colnames(data1) %in% c("treatment_outcome_first_course", "as.factor.myRsp.", "radiation_treatment_adjuvant", "targeted_molecular_therapy", "both_treatments", "progressive_disease")],
           interaction_col = "newTRT")

# Stable disease vs else (WITHOUT interaction terms)
fit_models(data2, ndpost = 5000, ntree = 1, seed = 42, 
           "stable_disease", 
           colnames(data2)[!colnames(data2) %in% c("treatment_outcome_first_course", "as.factor.myRsp.", "radiation_treatment_adjuvant", "targeted_molecular_therapy", "both_treatments", "progressive_disease", "stable_disease")],
           interaction_col = FALSE)

# Stable disease vs else (WITH interaction terms)
fit_models(data2, ndpost = 5000, ntree = 1, seed = 42, 
           "stable_disease", 
           colnames(data2)[!colnames(data2) %in% c("treatment_outcome_first_course", "as.factor.myRsp.", "radiation_treatment_adjuvant", "targeted_molecular_therapy", "both_treatments", "progressive_disease", "stable_disease")],
           interaction_col = "newTRT")

# Partial remission vs else (WITHOUT interaction terms)
fit_models(data3, ndpost = 5000, ntree = 1, seed = 42, 
           "partial_remission", colnames(data3)[!colnames(data3) %in% c("treatment_outcome_first_course", "as.factor.myRsp.", "radiation_treatment_adjuvant", "targeted_molecular_therapy", "both_treatments", "progressive_disease", "stable_disease", "partial_remission")],
           interaction_col = FALSE)

# Partial remission vs else (WITH interaction terms)
fit_models(data3, ndpost = 5000, ntree = 1, seed = 42, 
           "partial_remission", colnames(data3)[!colnames(data3) %in% c("treatment_outcome_first_course", "as.factor.myRsp.", "radiation_treatment_adjuvant", "targeted_molecular_therapy", "both_treatments", "progressive_disease", "stable_disease", "partial_remission")],
           interaction_col = "newTRT")

## -----------------------------------------------------------------------------


## Train/test of BART models ---------------------------------------------------
## Progressive vs else
data1 = read.csv("data1.csv")[,-(1:2)]
colnames(data1)
mean(data1$progressive_disease)
data_tree = data1[,!(colnames(data1)%in%c("treatment_outcome_first_course", "as.factor.myRsp.",
                                          "radiation_treatment_adjuvant", "targeted_molecular_therapy", "both_treatments",
                                          "progressive_disease"))]
# Fit the BART
set.seed(1)
n_train = floor(nrow(data_tree)*0.8)
bart_model = pbart(x.train=data_tree[1:n_train,], y.train=data1$progressive_disease[1:n_train],
                   ndpost=5000, ntree=1, base=0.95, power=0.5)
var_stable = bart_model$varcount.mean
sort(var_stable, decreasing = T)[1:5]
barplot(var_stable, main="Progressive vs. else")

# Posterior Predictive Distribution
y_pred = predict(bart_model, newdata=data_tree[(n_train+1):nrow(data_tree), ])
# yhat.test: Posterior predictive samples for the test data

y_pred_mean <- apply(y_pred$yhat.test, 2, mean)
y_pred_lower <- apply(y_pred$yhat.test, 2, quantile, probs = 0.025)
y_pred_upper <- apply(y_pred$yhat.test, 2, quantile, probs = 0.975)


prediction_intervals <- data.frame(
  y_pred_mean = y_pred_mean,
  y_pred_lower = y_pred_lower,
  y_pred_upper = y_pred_upper
)


# Create a data frame for plotting
plot_data <- data.frame(
  index = (n_train+1):nrow(data1),
  y_true = data1$progressive_disease[(n_train+1):nrow(data1)],
  y_pred_mean = y_pred_mean,
  lower = y_pred_lower,
  upper = y_pred_upper
)

ggplot(plot_data, aes(x = index, y = y_true)) +
  geom_point() +
  geom_line(aes(y = pnorm(y_pred_mean)), color = "blue") +
  geom_ribbon(aes(ymin = pnorm(lower), ymax = pnorm(upper)), alpha = 0.2) +
  labs(title = "BART Model Predictions with Uncertainty",
       x = "Index",
       y = "Value") +
  theme_minimal()


# set.seed(1)
# bart_model_sparse = pbart(x.train=data_tree, y.train=data1$progressive_disease,
#                     ndpost=5000, ntree=1, base=0.95, sparse=TRUE)
# sort(bart_model_sparse$varcount.mean, decreasing = T)[1:5]
# barplot(bart_model_sparse$varcount.mean)
# names(sort(bart_model_sparse$varcount.mean, decreasing = T)[1:5])


## Add interactions
data_tree_int = create_interaction_matrix(data_tree, "newTRT")
set.seed(1)
bart_model_int = pbart(x.train=data_tree_int, y.train=data1$progressive_disease,
                       ndpost=5000, ntree=1, base=0.9)
var_stable_int = bart_model_int$varcount.mean
barplot(var_stable_int, main="Progressive vs. else (with interactions)")
sort(var_stable_int, decreasing = T)[1:5]

## -----------------------------------------------------------------------------


## Tuning power and base parameter - 5-fold validation -------------------------
# Function for creating folds
create_folds = function(y, k){
  fold_indices = sample(rep(1:k, length.out = length(y)))
  folds = lapply(1:k, function(i) which(fold_indices == i))
  return(folds)
}

# Function for k-fold cross-validation
cv_bart2 = function(data, response_col, exclude_col, power_values, base_values, folds,
                    ndpost=5000, ntree=1){
  response = data[[response_col]]
  data_tree = data[,!(colnames(data)%in%exclude_col)]
  
  # Perform k-fold cross-validation
  k_folds = length(folds)
  cv_results = array(dim=c(length(power_values), length(base_values), k_folds),
                     dimnames=list(power_values, base_values, paste0("RMSE", 1:k_folds)))
  
  for(i in 1:length(power_values)){
    power = power_values[i]
    for(j in 1:length(base_values)){
      base = base_values[j]
      cat("\nTrying power", power,"and base", base,
          "---------------------------------------------------------------------\n")
      for(fold in 1:k_folds){
        # Split data into training and validation sets
        train_indices = sort(unlist(folds[-fold]))
        test_indices = folds[[fold]]
        
        x_train = data_tree[train_indices, ]
        y_train = response[train_indices]
        x_test = data_tree[test_indices, ]
        y_test = response[test_indices]
        
        # Fit the BART model with the current power value
        fit = pbart(x.train=x_train, y.train=y_train,
                    ndpost=ndpost, ntree=ntree, base=base, power=power)
        
        # Generate predictions and calculate RMSE for the current fold
        y_pred = predict(fit, newdata=x_test)$yhat.test
        y_pred_mean = apply(y_pred, 2, mean)
        rmse = sqrt(mean((y_test - y_pred_mean)^2))
        cv_results[i,j,fold] = rmse
      }
    }
  }
  
  return(cv_results)
}

power_values = seq(0.1, 2.5, by=0.2) # beta
base_values = seq(0.1, 2, by=0.2) # alpha
k_folds = 5
set.seed(1)
folds = create_folds(data1$progressive_disease, k_folds)

# Cross-validation results of progressive_disease
data1 = read.csv("data1.csv")[,-(1:2)]
cv_results = cv_bart2(data=data1, response_col="progressive_disease", 
                      exclude_col=c("treatment_outcome_first_course", "as.factor.myRsp.",
                                    "radiation_treatment_adjuvant", "targeted_molecular_therapy", "both_treatments",
                                    "progressive_disease"),
                      power_values=power_values, base_values=base_values, folds=folds,
                      ndpost=5000, ntree=1)
cv_results_mean = apply(cv_results, 1:2, mean)
cv_results_mean

cv_results_mean_df = as.data.frame(as.table(cv_results_mean))
colnames(cv_results_mean_df) = c("X", "Y", "Z")
ggplot(cv_results_mean_df, aes(x=X, y=Y, fill=Z)) +
  geom_tile() +
  labs(title = "Cross-Validation for 'No Progression vs. Else' - Bayesian CART",
       x = "Power",
       y = "Base",
       # z = "RMSE",
       fill ="RMSE") +
  scale_fill_viridis_c(direction = +1)#, option="B"


# Cross-validation results of stable_disease
data2 = read.csv("data2.csv")[,-(1:2)]
set.seed(1)
folds = create_folds(data2$stable_disease, k_folds)
cv_results = cv_bart2(data=data2, response_col="stable_disease",
                      exclude_col=c("treatment_outcome_first_course", "as.factor.myRsp.",
                                    "radiation_treatment_adjuvant", "targeted_molecular_therapy", "both_treatments",
                                    "progressive_disease", "stable_disease"),
                      power_values=power_values, base_values=base_values, folds=folds,
                      ndpost=5000, ntree=1)
cv_results
# Summarize the cross-validation results
cv_results_mean = apply(cv_results, 1, mean)
plot(cv_results_mean, type="b", pch=16, xaxt="n",
     ylim=c(min(cv_results_mean)-0.01, max(cv_results_mean)+0.01),
     ylab="RMSE", xlab="Power", main="Cross-Validation Results for Bayesian CART Power Parameter")
axis(side=1, at=1:length(power_values), labels=power_values)
grid(nx=NULL, ny=NULL, col="gray", lty=3,  lwd=1)


## -----------------------------------------------------------------------------


