#' lmeBioCov
#'
#' This function fits linear mixed-effects models to test the influence of a
#' set of covariates in the levels of the dependent variable (biomarkers). It takes
#' as input a dataframe ("df") containing the data to be analyzed, a vector with
#' the names of the biomarkers ("bio_names") to be tested, a vector of individual
#' covariates ("ind_cov") in which we will test the interaction, the covariate
#' ("cov") that we are adding to the model besides the individual covariates,
#' the name of the ID or patient column in the top-table ("patientCol"), the
#' formula of the fixed effects of the model ("fixedFormula"), the formula of the
#' random effects of the model ("randomFormula"), and the algorithm used for
#' estimating the model parameters (REML). The function returns a dataframe
#' with the results of the different linear models for each biomarker, including
#' the coefficients, standard errors, t-values, p-values, R-squared, AIC, and
#' the number of outliers.
#'
#' @param df A given dataframe
#' @param bio_names Vector with the names of the biomarkers to test.
#' @param ind_cov Individual covariates for which we want to test the interaction
#' with the biomarkers. They will be included as fixed effects in all models.
#' @param covariate The additional covariate that we are including in the model.
#' @param patientCol Name of the ID or patient column in the top-table.
#' @param fixedFormula Formula of the Fixed Effects of the Model.
#' @param randomFormula Formula of the Random Effects of the Model.
#' @param REML Algorithm used for estimating the model parameters. It is recommended
#' to use FALSE (it will use ML) if we are comparing models.
#' @return A dataframe with the results of the different linear models.
#' @import lmerTest rsq progress
#' @importFrom dplyr bind_rows
#' @export
#'

lmeBioCov <- function(df, bio_names, ind_cov, covariate, patientCol, fixedFormula, randomFormula, REML){
  # First, we get the indexes of the biomarkers, patient and individual covariates:
  vars <- names(df)
  bio_id <- which(vars %in% bio_names)
  patient_id <- which(vars == patientCol)
  cov_id <- which(vars %in% ind_cov)
  var_id <- which(vars %in% covariate)

  # And we create a dataframe to store the data:
  lme_df <- data.frame()

  # We start a counter:
  j <- 0

  # We create a progress bar to show the time remaining for every model:
  pb <- progress_bar$new(total = length(bio_id), format = "[:bar] :percent :eta")

  # And we loop through all of the biomarkers:
  for (i in bio_id){
    # First, we add to the counter:
    j <- j + 1

    # We create a subset with the target columns: biomarker, patient and covariates:
    df_bio <- df[, c(patient_id, i, cov_id, var_id)]

    # We create a dataframe without NA to calculate the R2 coefficient:
    df_clean <- na.omit(df_bio)

    # We change the name of the biomarker column to "biomarker" to have the same
    # name as in the model formula:
    colnames(df_bio)[2] <- "biomarker"

    # Now, we generate the model:
    fit <- lmer(formula = paste(fixedFormula, "+", randomFormula), data = df_bio, REML = REML)
    # And we extract the summary of the model:
    sum_table <- data.frame(coef(summary(fit)))

    # We change the name of the p-value column:
    colnames(sum_table)[5] <- "p-value"

    # We compute the R2 coefficient of the model and its components using rsq():
    R2 <- rsq(fit)

    # And then we format our summary table with the model information:
    sum_table$biomarker <- names(df)[i]
    sum_table$coefficient <- rownames(sum_table)
    sum_table$R2_model <- R2$model
    sum_table$R2_fixed_effects <- R2$fixed
    sum_table$R2_random_effects <- R2$random

    # And we also extract the AIC and BIC parameters from the summary:
    sum_table$AIC <- summary(fit)$AICtab[1] #AIC
    sum_table$BIC <- summary(fit)$AICtab[2] #BIC

    # Finally, we add the rows generated for each protein:
    lme_df <- bind_rows(lme_df, sum_table)

    # And we add a tick to the progress bar:
    pb$tick()
  }

  # Finally, we return the dataframe with the information:
  return(lme_df)

}
