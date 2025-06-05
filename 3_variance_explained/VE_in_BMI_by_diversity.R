#Mariam Miari
#250604
# variance explained by diversity metrics in BMI

rm(list = ls())
library(rio)
library(BiocParallel) 
library(glmnet)
library(caret)
library(tidyverse)
library(readxl)
library(here)
library(vegan)
library(openxlsx) # for exporting results

set.seed(1)
load("../data/processed/MGS2.0/data_ds_nontrans_EUR_240320.RData")
data_ds_eur <- data_ds_eur[complete.cases(data_ds_eur$UNRELATED_1_DEGREE),]


mgs_st <- which(colnames(data_ds_eur) == "shannon_ds") +1
mgs_nontrans <- data_ds_eur[,c(mgs_st:ncol(data_ds_eur))]

#calc inverse simpson
inverse_simpson_values <- diversity(mgs_nontrans, index = "invsimpson", MARGIN = 1)
inverse_simpson_values <- as.data.frame(cbind(data_ds_eur$Data_sample_code, inverse_simpson_values))
colnames(inverse_simpson_values)[1] <- "Data_sample_code"

merged_df <- merge(inverse_simpson_values, data_ds_eur, by = "Data_sample_code")

merged_df$sex <- ifelse(merged_df$sex == 1,0,1)

#calc richness
richness_ds <- as.data.frame(rowSums(mgs_nontrans > 0))
colnames(richness_ds)[1] <- "richness_ds"
rownames(richness_ds) <- merged_df$Data_sample_code

merged_df <- merge(richness_ds, merged_df, by.x = "row.names", by.y = "Data_sample_code")

# get variables
tax <- merged_df[,c("richness_ds","shannon_ds","inverse_simpson_values",'age', "sex")]
rownames(tax) <- merged_df$Data_sample_code
tax$inverse_simpson_values <- as.numeric(tax$inverse_simpson_values)

pheno <- merged_df$BMI

######################################################## all individuals
# function to compute R² for a given diversity metric
var.fun <- function(pheno, tax, metric) {
  
  data <- data.frame(pheno = pheno, tax)
  
  # rename the selected diversity metric dynamically
  colnames(data)[which(colnames(data) == metric)] <- "diversity_metric"
  
  # create 10 cross-validation folds
  cv_folds <- createFolds(data$pheno, k = 10)
  
  cv <- lapply(cv_folds, function(cv_fold) {
    
    train <- data[-cv_fold, ]  # 9 folds for training
    test <- data[cv_fold, ]    # 1 fold for testing
    
    # full model (with diversity metric)
    fit_full <- glm(pheno ~ diversity_metric + age + sex, data = train, family = gaussian())
    
    # null model (without diversity metric)
    fit_null <- glm(pheno ~ age + sex, data = train, family = gaussian())
    
    pred_full_test <- predict(fit_full, test, type = "response")
    pred_null_test <- predict(fit_null, test, type = "response")
    
    # calculate R² using test set
    r.squared.full.test <- cor(pred_full_test, test$pheno)^2
    r.squared.null.test <- cor(pred_null_test, test$pheno)^2
    
    # R² difference (contribution of diversity metric)
    r.squared.diff.test <- r.squared.full.test - r.squared.null.test
    
    data.frame(r.squared.full.test, r.squared.null.test, r.squared.diff.test)
  })
  
  cv <- do.call(rbind, cv)
  
  # return mean R² across all folds
  data.frame(
    metric = metric,  # store the diversity metric name
    r.squared.full.test = mean(cv$r.squared.full.test, na.rm = TRUE), 
    r.squared.null.test = mean(cv$r.squared.null.test, na.rm = TRUE),
    r.squared.diff.test = mean(cv$r.squared.diff.test, na.rm = TRUE)
  )
}

diversity_metrics <- c("shannon_ds", "inverse_simpson_values", "richness_ds")

# run the function on non-bootstrapped results for each diversity metric
var_exp_est <- do.call(rbind, lapply(diversity_metrics, function(metric) {
  var.fun(pheno, tax, metric)
}))

# bootstrapping
n_boot <- 500  
n <- length(pheno)  
var_boot_ls <- vector("list", length = n_boot)

set.seed(123)

for(i in 1:n_boot){
  
  boot_ix <- sample.int(n, replace = TRUE)  
  
  pheno_boot <- pheno[boot_ix]   
  tax_boot <- tax[boot_ix, ]
  
  # run function for all diversity metrics
  var_boot_ls[[i]] <- do.call(rbind, lapply(diversity_metrics, function(metric) {
    var.fun(pheno_boot, tax_boot, metric)
  }))
}

# combine bootstrapped results into one dataframe
combined_df <- do.call(rbind, var_boot_ls)

# standard deviation of bootstrapped R² values
confidence_level <- 0.975
critical_value <- qnorm(confidence_level)

var_exp_var <- combined_df %>%
  group_by(metric) %>%
  summarise(
    variance = var(r.squared.diff.test, na.rm = TRUE),
    stdev = sqrt(variance),
    .groups = "drop"
  )

# merge bootstrapped variance with non-bootstrapped results
var_exp_all <- merge(var_exp_est, var_exp_var, by = "metric")

# calculate confidence intervals
var_exp_all <- var_exp_all %>% 
  mutate(
    upper_CI = r.squared.diff.test + critical_value * stdev,
    lower_CI = r.squared.diff.test - critical_value * stdev,
    analysis = "all"
  )

print(var_exp_all)

# export results
write.xlsx(var_exp_all, file = "results/VE_by_diversity_metrics_with_CI.xlsx")

