---
title: "alpha_diversity_regression"
output: html_document
date: "2025-06-05"
author: Mariam Miari
description: 1. Association of alpha diversity metrics with BMI as outcome
              2. Association of BMI with alpha diversity metrics as outcome
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

#load libraries
```{r}
rm(list = ls())
library(tidyverse)
library(ggplot2)
library(here)
library(readxl)
library(vegan)
library(gridExtra)  
library(scales)    
library(compositions)
library(reshape2)
library(rstatix)
library(rio)
```

#load data
```{r}
data_trans <- read_excel("../data/processed/MGS2.0/data_nds_filt_9993_trans_EUR_PRS_240624.xlsx") # this df has shannon_nds and shannon_ds

# downsized non transformed
load("../data/processed/MGS2.0/data_ds_nontrans_EUR_240320.RData") #has downsized shannon
data_ds_eur <- data_ds_eur[complete.cases(data_ds_eur$UNRELATED_1_DEGREE),]

```

##preprocessing
```{r}
mgs_st <- which(colnames(data_ds_eur) == "shannon_ds") +1
mgs_nontrans <- data_ds_eur[,c(mgs_st:ncol(data_ds_eur))]
richness_ds <- as.data.frame(rowSums(mgs_nontrans > 0))
colnames(richness_ds)[1] <- "richness_ds"
richness_df <- cbind(data_ds_eur$Data_sample_code, richness_ds)
colnames(richness_df)[1] <- "Data_sample_code"
```

```{r}
inverse_simpson_values <- diversity(mgs_nontrans, index = "invsimpson", MARGIN = 1)
inverse_simpson_df <- as.data.frame(cbind(data_ds_eur$Data_sample_code, inverse_simpson_values))
colnames(inverse_simpson_df)[1] <- "Data_sample_code"
inverse_simpson_df$inverse_simpson_values <- as.numeric(inverse_simpson_df$inverse_simpson_values)
```

```{r}
Data_sample_code <- data_ds_eur[,c("Data_sample_code")]
div_df <- merge(richness_df, inverse_simpson_df, by = "Data_sample_code")
```


```{r}
data <- merge(div_df, data_ds_eur, by = "Data_sample_code")

data <- data %>%
  mutate(physical_activity_f = factor(physical_activity)) %>%
  relocate(physical_activity_f, .after = physical_activity)

#factor smoking
data <- data %>%
  mutate(smoking_f = factor(smoking)) %>%
  relocate(smoking_f, .after = smoking)

data <- data %>% 
  mutate(male = sex == 1) %>%
  relocate(male, .after = sex)
```


##regression
###BMI as outcome association with diversity metrics
```{r}
# alpha diversity as predictor and BMI as outcome
run_glm_div_as_exposure <- function(diversity_metric, data, sex = "all", interaction = FALSE) {
  
  confounders <- c("age", "country_birth", "Cohort", "physical_activity_f",
                   "smoking_f", "metformin", "diab_diag", "Fiber")
  
  if (sex == "all") {
    confounders <- c("male", confounders)
  }
  
  if (interaction) {
    formula <- as.formula(paste("BMI ~", diversity_metric, "* male +", paste(confounders, collapse = " + ")))
  } else {
    formula <- as.formula(paste("BMI ~", diversity_metric, "+", paste(confounders, collapse = " + ")))
  }
  
  #  assign model_name
  if (sex == "male") {
    data <- subset(data, male == TRUE)
    model_name <- "sex_stratification-male"
  } else if (sex == "female") {
    data <- subset(data, male == FALSE)
    model_name <- "sex_stratification-female"
  } else {
    model_name <- "all"
  }
  
  if (interaction) {
    model_name <- "interaction"
  }
  
  model <- glm(formula, data = data, family = gaussian())
  
  coef_table <- summary(model)$coefficients
  term <- if (interaction) paste0(diversity_metric, ":maleTRUE") else diversity_metric
  
  if (!(term %in% rownames(coef_table))) return(NULL)
  
  coef_row <- coef_table[term, ]
  conf_int <- confint(model, term)
  
  results <- data.frame(
    exposure = diversity_metric,
    model_name = model_name,
    n_obs = model$df.null + 1,
    est = coef_row[1],
    se = coef_row[2],
    ci_low = conf_int[1],
    ci_high = conf_int[2],
    p = coef_row[4]
  )
  
  return(results)
}

#diversity metrics as predictors of BMI
diversity_metrics <- c("shannon_ds", "richness_ds", "inverse_simpson_values")
results_list_bmi <- list()

for (metric in diversity_metrics) {
  results_list_bmi[[length(results_list_bmi) + 1]] <- run_glm_div_as_exposure(metric, data, sex = "all")
  results_list_bmi[[length(results_list_bmi) + 1]] <- run_glm_div_as_exposure(metric, data, sex = "male")
  results_list_bmi[[length(results_list_bmi) + 1]] <- run_glm_div_as_exposure(metric, data, sex = "female")
  results_list_bmi[[length(results_list_bmi) + 1]] <- run_glm_div_as_exposure(metric, data, sex = "all", interaction = TRUE)
}

results_bmi_df <- do.call(rbind, results_list_bmi)

export(results_bmi_df, "v3/results_regression_diversity_BMI_sex_strat_int.xlsx")

```

###BMI as exposure association with diversity metrics (outcome)
```{r}
run_glm_bmi_as_exposure <- function(outcome_metric, data, sex = "all", interaction = FALSE) {
  
  confounders <- c("age", "country_birth", "Cohort", "physical_activity_f",
                   "smoking_f", "metformin", "diab_diag", "Fiber")
  
  if (sex == "all") {
    confounders <- c("male", confounders)
  }
  
  if (interaction) {
    formula <- as.formula(paste(outcome_metric, "~ BMI * male +", paste(confounders, collapse = " + ")))
  } else {
    formula <- as.formula(paste(outcome_metric, "~ BMI +", paste(confounders, collapse = " + ")))
  }
  
  # data and assign model_name
  if (sex == "male") {
    data <- subset(data, male == TRUE)
    model_name <- "sex_stratification-male"
  } else if (sex == "female") {
    data <- subset(data, male == FALSE)
    model_name <- "sex_stratification-female"
  } else {
    model_name <- "all"
  }
  
  if (interaction) {
    model_name <- "interaction"
  }
  
  model <- glm(formula, data = data, family = gaussian())
  
  coef_table <- summary(model)$coefficients
  term <- if (interaction) "BMI:maleTRUE" else "BMI"
  
  if (!(term %in% rownames(coef_table))) return(NULL)
  
  coef_row <- coef_table[term, ]
  conf_int <- confint(model, term)
  
  results <- data.frame(
    exposure = "BMI",
    outcome = outcome_metric,
    model_name = model_name,
    n_obs = model$df.null + 1,
    est = coef_row[1],
    se = coef_row[2],
    ci_low = conf_int[1],
    ci_high = conf_int[2],
    p = coef_row[4]
  )
  
  return(results)
}

# BMI as predictor of diversity metrics
diversity_metrics <- c("shannon_ds", "richness_ds", "inverse_simpson_values")
results_list <- list()

for (metric in diversity_metrics) {
  results_list[[length(results_list) + 1]] <- run_glm_bmi_as_exposure(metric, data, sex = "all")
  results_list[[length(results_list) + 1]] <- run_glm_bmi_as_exposure(metric, data, sex = "male")
  results_list[[length(results_list) + 1]] <- run_glm_bmi_as_exposure(metric, data, sex = "female")
  results_list[[length(results_list) + 1]] <- run_glm_bmi_as_exposure(metric, data, sex = "all", interaction = TRUE)
}

results_bmi_as_exposure_df <- do.call(rbind, results_list)

export(results_bmi_as_exposure_df, "v3/results_regression_BMI_diversity_sex_strat_int.xlsx")

```

