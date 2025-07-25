# Marlena Maziarz
# June 8, 2022, edited August 18, 23, 30, September 29, 2022
# for the microbiome BMI paper

# functions used by the main script 'bmi_paper_model_fitting.Rmd'

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function definitions
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# helper functions
get_model_formula_helper <- function(outcome,                      # string
                                     mgs,                          # string
                                     base_model,                   # T/F
                                     base_model_vars_to_adjust_for,
                                     full_model_vars_to_adjust_for,
                                     mgs_sex_interaction,
                                     sex_var){
  
  if(base_model){
    vars_to_adjust_for <- c(mgs, base_model_vars_to_adjust_for)
  }else{
    vars_to_adjust_for <- c(mgs, full_model_vars_to_adjust_for)
  }
  
  model_formula <- paste0(outcome, ' ~ ', paste(vars_to_adjust_for, collapse = ' + '))
  
  # add the interaction term if the model is specified with an interaction
  if(mgs_sex_interaction){
    model_formula <- paste0(model_formula, ' + ', mgs, ':', sex_var)
  }
  
  return(as.formula(model_formula))
}



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# the main function
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# NOTE: the description for this function has not yet been updated after changes to the function
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# input:
# d = name of the dataset with samples in rows, MGSs and other covariates in columns
#     this is the only one with no default
#
# outcome = c('BMI'), # outcome of interest
#
#
# output:
# a data frame that looks like this:

# > head(final_results)
#   model           mgs         est        se          p        p_adj
# 1 m_base_BMI     HG3A.0001 -0.22760398 0.2642608 0.38909900    NA
# 2 m_base_BMI HG3A.0001:sex  0.26238741 0.1637971 0.10920427    NA
# 3 m_base_BMI     HG3A.0002 -1.19539864 0.5481074 0.02920701    NA
# 4 m_base_BMI HG3A.0002:sex  0.59008442 0.3369835 0.07995945    NA
# 5 m_base_BMI     HG3A.0003  0.06488077 0.4336367 0.88106700    NA
# 6 m_base_BMI HG3A.0003:sex -0.03846190 0.2687106 0.88618608    NA
# ...
# currently it outputs results for 6 models, base and full for the three outcomes.
# as it runs, it prints each of the 6 models as well as the estimates for the models for the first MGS.

get_results_for_bmi_paper <- function(d,       # dataset must be specified
                                      outcome, # the outcome must be specified 
                                      model_name = NULL,
                                      model_strata = 'all', # 'all', 'males', 'females'
                                      suppress_printing = F,
                                      base_model_vars_to_adjust_for = c('sex', 'age', 'country_birth', 'Cohort'), # base model must be provided
                                      full_model_vars_to_adjust_for = c('physical_activity_f', 'smoking_f', 'metformin', 'diab_diag', 'Fiber'), # can be NA
                                      glm_family = 'gaussian', #with glmrob() argument method = "Mqle" is added
                                      mgs_sex_interaction = F,
                                      sex_var = 'male'){ # for interactions with sex and stratification based on sex
  
  d <- d %>%
    droplevels() # make sure to get rid of all the 'plate', etc. levels that are not in the dataset 
  
  # get the names of all MGS's in the dataset
  all_mgs <- names(d)[grep('hMGS', names(d))]
  n_mgs = length(all_mgs)
  
  # 6-10 depending on the bmi_adj_for_waist and bmi_adj_for_WHR flags
  nrows <- n_mgs * (length(outcome)) * 
    (1 + ifelse(is.na(full_model_vars_to_adjust_for), 0, 1)) * 
    (1 + mgs_sex_interaction)
  
  results_ix <- 1:(1 + mgs_sex_interaction) # 1:1 or 1:2, the initial index for the results in the results matrix
  index_increment <- 1 + mgs_sex_interaction  # 1 or 2
  
  
  # make a 'dummmy' matrix for the results
  results_columns <- c('model_name', # model name
                       'mgs',   # name of the MGS, or MGS:male, etc
                       'n_obs', # the number of observations used in the model fitting (from the glm object)
                       'est', 'se', 'ci_low', 'ci_high', 'p', 'p_fdr_per_outcome')
  if(glm_family == 'binomial'){
    # make sure the user knows that the results are not exponentiated.
    cat('\n\n***** NOTE: The est in the results matrix are log(odds ratios), use exp(est) to get odds ratios. *****\n\n')
  }
  
  results <- as.data.frame(matrix(NA, nrow = nrows, ncol = length(results_columns)))
  colnames(results) <- results_columns
  
  for(base_model in c(T, F)){
    
    if(is.null(model_name)){
      # default name
      curr_model_name <- ifelse(base_model, paste0('m_base_', outcome), paste0('m_full_', outcome))
    }else{
      curr_model_name <- ifelse(base_model, paste0('base_', model_name), paste0('full_', model_name))
    }
    
    
    for (i in 1:n_mgs){
      if(i == 1 & suppress_printing == F){
        print_curr_model <- T # print the model specification once for each type of model
      }else{
        print_curr_model <- F
      }
      curr_model_formula <- get_model_formula_helper(outcome = outcome, 
                                                     mgs = all_mgs[i], 
                                                     base_model = base_model, # T/F
                                                     base_model_vars_to_adjust_for = base_model_vars_to_adjust_for,
                                                     full_model_vars_to_adjust_for = full_model_vars_to_adjust_for,
                                                     mgs_sex_interaction = mgs_sex_interaction,
                                                     sex_var = sex_var)
      #notice glm instead of glm, this is glm robust regression to handle leveraging in clr transfrmation
      if(model_strata == 'all'){
        curr_model <- glm(curr_model_formula, family = glm_family, data = d)  #with glmrob() argument method = "Mqle" is added 
      }else if(model_strata == 'males'){
        curr_model <- glm(curr_model_formula, family = glm_family, data = d, subset = get(sex_var) == T) #with glmrob() argument method = "Mqle" is added
      }else if(model_strata == 'females'){
        curr_model <- glm(curr_model_formula, family = glm_family, data = d, subset = get(sex_var) == F) #with glmrob() argument method = "Mqle" is added
      }else{
        print('No such stratification option - your choices are: all, males, females. Exiting...')
        break
        
      }
      
      # print the model for the first MGS, just to verify that we're fitting the models we wanted
      # only print the model summary if plate is not a covariate
      if(print_curr_model & length(grep('plate', curr_model_formula)) == 0){
        cat(paste('~~~~~~~~~~~~~~~~~ start of model', curr_model_name, '~~~~~~~~~~~~~~~~~~~\n\n'))
        cat(paste('Model formula for', curr_model_name, ':\n'))
        print(curr_model_formula, showEnv = F, q = F)
        cat(' \nModel estimates for the first MGS:\n') 
        print(summary(curr_model)$coeff, q = F) 
        cat(paste('~~~~~~~~~~~~~~~~~~ end of model', curr_model_name, '~~~~~~~~~~~~~~~~~~~~ \n\n'))
      }
      
      # print an error if any of the estimates are missing (NAs)
      if(sum(is.na(summary(curr_model)$coeff)) > 0){
        print('ERROR: there are NAs in the model. Skipping this model and moving on to the next.', q = F)
        print('Info about the model that generated the NAs:', q = F)
        print(curr_model_name, q = F)
        print(curr_model_formula, q = F)
        print(summary(curr_model)$coeff, q = F)
        break
      }
      
      est_of_interest_row_ix <- grep(all_mgs[i], names(curr_model$coeff)) # this will be wrong if NA's in the model
      if(glm_family == 'gaussian'){
        est_of_interest <- summary(curr_model)$coeff[est_of_interest_row_ix, c('Estimate', 'Std. Error', 'Pr(>|t|)')] #with glmrob() Pr(>|z|)
      }else if(glm_family == 'binomial'){
        est_of_interest <- summary(curr_model)$coeff[est_of_interest_row_ix, c('Estimate', 'Std. Error', 'Pr(>|z|)')]
      }else{
        print(paste0('ERROR: there is no support for glm_family = ', glm_family, '. Returning NULL.'), q = F)
        return(NULL)
      }
      
      # save the results in the 'final' results matrix
      results[results_ix, 'model_name'] <- curr_model_name
      results[results_ix, 'mgs'] <- names(curr_model$coeff)[est_of_interest_row_ix] # this works even if est is NA
      results[results_ix, 'n_obs'] <- curr_model$df.null +1 # with glmrob() length(curr_model$y)
      results[results_ix, c('est', 'se', 'p')] <- est_of_interest
      
      # update the index for the next iteration
      results_ix <- results_ix + index_increment
      
      # cleanup (not needed if all goes well, but will reveal a problem if something goes wrong,
      #          if, for example, the model doesn't fit for a given MGS)
      rm(curr_model_formula, curr_model, est_of_interest_row_ix, est_of_interest)
    }
    
  }
  
  # calculate the confidence intervals
  results[, c('ci_low', 'ci_high')] <- with(results, est + qnorm(0.975) * se %*% t(c(-1, 1)))
  
  
  return(results)
}

# if you run into any problems, try to find the error using the debug() function
# debug('get_results_for_bmi_paper')

