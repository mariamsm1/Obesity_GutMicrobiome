######### Mariam Miari
######### 250604
#########LASSO regression to get variance explained by species in BMI
######### The same script is run for sex-stratified analysis. Subset the data accordingly.

rm(list=ls())
library(rio)
library(BiocParallel) 
library(glmnet)
library(caret)
library(tidyverse)
library(ggplot2)
library(readxl)
library(here)


set.seed(1)

##################### lasso regression 
## 
data <- readxl::read_xlsx("../data/processed/MGS2.0/data_nds_filt_9993_trans_EUR_PRS_240624.xlsx") 

  
data <- data %>% 
  mutate(male = sex == 1)


rownames(data) <- data$Data_sample_code
BMI <- data$BMI


## all phenotypes
pheno <- as.data.frame(BMI)


## Get the MGSs
start <- which(colnames(data) == 'shannon_nds')+1
tax <- data[start:(ncol(data)-1)] 


var.fun <- function(pheno, tax) {
  
  res <- bplapply(colnames(pheno), function(x) {
    
    # create data for complete cases
    
    data <- data.frame(pheno = pheno[, x], tax)
    pheno <- data$pheno
    tax <- data[, which(!colnames(data) == "pheno")] #####
    
    tryCatch({
      
      # create 10 folds
      
      folds <- createFolds(pheno, k = 10)
      
      nested.cv <- lapply(folds, function(fold) {
        
        # create test data
        
        pheno_test <- pheno[fold] # PHENOTYPE
        tax_test <- tax[fold, ] # MGS
        
        # create cross-validation data
        
        pheno_cv <- pheno[-fold]
        tax_cv <- tax[-fold, ]
        
        # create 10 cross-validation folds
        
        cv_fold <- createFolds(pheno_cv, k = 10)[[1]]
        
        # create validation data
        
        pheno_validate <- pheno_cv[cv_fold]
        tax_validate <- tax_cv[cv_fold, ]
        
        # create training data
        
        pheno_train <- pheno_cv[-cv_fold]
        tax_train <- tax_cv[-cv_fold, ]
        
        # fit lasso regression model
        
        fit <- glmnet(x = as.matrix(tax_train), y = pheno_train, alpha = 1)
        
        # determine model with minimum lambda
        
        pred <- predict(fit, newx = as.matrix(tax_validate), s = fit$lambda)
        mse <- colMeans((pred - pheno_validate)^2)
        lambda <- fit$lambda[which.min(mse)]
        mse.train <- min(mse)
        
        # test model with minimum lambda on test fold
        
        pred <- as.numeric(predict(fit, newx = as.matrix(tax_test), s = lambda))
        r.squared <- cor(pred, pheno_test)^2
        mse.test <- mean((pred - pheno_test)^2)
        data.frame(r.squared = r.squared, mse.train = mse.train, mse.test = mse.test)
        
      })
      nested.cv <- do.call(rbind, nested.cv)
      
      
      data.frame(phenotype = x, r.squared = mean(nested.cv$r.squared), mse.train = mean(nested.cv$mse.train), mse.test = mean(nested.cv$mse.test), n = nrow(tax), error = NA)
    }, error = function(e) {
      
      # return empty result if error
      
      data.frame(phenotype = NA, r.squared = NA, mse.train = NA, mse.test = NA, n = NA, error = paste("Error:", e$message))
      
    })
    
  }, BPPARAM = MulticoreParam(workers = 16))
  
  do.call(rbind, res)
  
}

save(pheno, tax, file = here::here('data.rda'))

# import data

load(here::here("data.rda"))



variance.explained <- var.fun(pheno,tax)

# order

variance.explained <- variance.explained[order(-abs(variance.explained$r.squared)), ]

# export data

export(variance.explained[, c("phenotype", "r.squared", "mse.train", "mse.test", "n")],
       here::here("results","variance_explained_MGS_resid_lasso.tsv"))


########################### bootstraping for generating CIs
bootstrap_var_fun <- function(pheno, tax, n_bootstrap = 500, conf_level = 0.95) {
  
  # calculate original estimates on the full dataset
  original_est <- var.fun(pheno, tax)
  
  # bootstrap resampling
  bootstrap_results <- replicate(n_bootstrap, {
    
    # sample with replacement from indices of the data
    sample_indices <- sample(1:nrow(pheno), size = nrow(pheno), replace = TRUE)
    pheno_sample <- pheno[sample_indices, , drop = FALSE]
    tax_sample <- tax[sample_indices, , drop = FALSE]
    
    # run var.fun on the resampled data
    tryCatch({
      var.fun(pheno_sample, tax_sample)
    }, error = function(e) {
      # return NA data frame if error occurs
      data.frame(phenotype = NA, r.squared = NA, mse.train = NA, mse.test = NA, n = NA, error = paste("Error:", e$message))
    })
    
  }, simplify = FALSE)
  
  # combine results
  bootstrap_results <- do.call(rbind, bootstrap_results)
  
  bootstrap_results <- bootstrap_results[!is.na(bootstrap_results$r.squared), ]
  
  # check if bootstrap results are empty after filtering
  if (nrow(bootstrap_results) == 0) {
    stop("All bootstrap samples resulted in errors or NA values. Check your `var.fun` function or data for potential issues.")
  }
  
  # convert columns to numeric to ensure calculations work smoothly
  bootstrap_results$r.squared <- as.numeric(bootstrap_results$r.squared)
  bootstrap_results$mse.train <- as.numeric(bootstrap_results$mse.train)
  bootstrap_results$mse.test <- as.numeric(bootstrap_results$mse.test)
  
  # calculate variance and standard deviation of r.squared across bootstraps
  var_exp_var <- bootstrap_results %>%
    group_by(phenotype) %>%
    summarise(
      variance = var(r.squared, na.rm = TRUE),
      stdev = sqrt(variance)
    )
  
  # merge with the original estimate
  combined_results <- merge(var_exp_var, original_est, by = "phenotype")
  
  # calculate the critical value for the desired confidence level
  critical_value <- qnorm((1 + conf_level) / 2)
  
  # calculate confidence intervals for r.squared based on the original estimate
  combined_results <- combined_results %>%
    mutate(
      r.squared = r.squared,  # original r.squared from the non-bootstrapped model
      upper_CI = r.squared + critical_value * stdev,
      lower_CI = r.squared - critical_value * stdev
    ) %>%
    select(phenotype, r.squared, variance, stdev, lower_CI, upper_CI, mse.train, mse.test, n)
  
  return(combined_results)
}


bootstrap_results <- bootstrap_var_fun(pheno, tax)
export(bootstrap_results, here::here("results/VE_BMI_with_CI.xlsx"))

