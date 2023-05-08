#' lmeBioInt
#'
#' This function, named "lmeBioInt", is a part of a package designed to perform
#' linear mixed-effects modeling in R. The purpose of this function is to fit a
#' linear mixed-effects model to test the interaction between given biomarkers
#' and individual covariates (with the optional addition of other covariates in
#' the model) in a given data frame. The function takes a data frame as input
#' along with the names of the biomarkers and individual covariates to test,
#' the name of the patient column, the dependent variable, the covariates
#' to include in the model (optional, NULL by default), the formulas for fixed
#' and random effects, and the method used for estimating the model parameters (REML).

#' The function first identifies the indexes of the biomarkers, patient, and
#' individual covariates and creates a data frame to store the results. Then,
#' the function loops through each biomarker and generates a linear model with
#' the provided information. Then, it extracts the summary table and creates a
#' dataframe with model information, including the predicted value of the model,
#' marginal and conditional residuals, and outliers. Finally, the function
#' calculates the AIC and returns a data frame with all the information. The
#' function uses the "lme" function from the "nlme" package to fit the model.
#'
#' @param df A given dataframe
#' @param bio_names Vector with the names of the biomarkers to test.
#' @param ind_cov Individual covariates for which we want to test the interaction
#' with the biomarkers. They will be included as fixed effects in all models.
#' @param patientCol Name of the ID or patient column in the top-table.
#' @param depVar Dependent variable, as clinical outcome or score.
#' @param covariates Vector with the covariates that we want to include in the model
#' (NULL by default to test interaction of individual covariates with biomarkers).
#' @param fixedFormula Formula of the Fixed Effects of the Model.
#' @param randomFormula Formula of the Random Effects of the Model.
#' @param REML Algorithm used for estimating the model parameters. It is recommended
#' to set it as FALSE (it will use ML) if we are comparing models.
#' @return A dataframe with the results of the different linear models.
#' @importFrom dplyr bind_rows
#' @import openxlsx lmerTest rsq progress
#' @rawNamespace import(kableExtra, except = group_rows)
#' @export
#'

lmeBioInt <- function(df, bio_names, ind_cov, patientCol, depVar, covariates = NULL, fixedFormula, randomFormula, REML){
  # First, we get the indexes of the biomarkers, patient and individual covariates:
  vars <- names(df)
  bio_id <- which(vars %in% bio_names)
  patient_id <- which(names(df) == patientCol)
  dep_id <- which(names(df) == depVar)
  cov_id <- which(vars %in% ind_cov)
  covariate_id <- which(names(df) == covariates)

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
    df_bio <- df[, c(patient_id, dep_id, i, cov_id, covariate_id)]

    # We adapt the fixed formula depending on if we are adding a covariate to the
    # model of if we are only testing the interaction:
    if (!is.null(covariates)) {
      fixedFormula <- paste(depVar, "~", paste0(ind_cov, collapse = "*"), "*", colnames(df_bio)[3], "+", covariates)
    } else {
      fixedFormula <- paste(depVar, "~", paste0(ind_cov, collapse = "*"), "*", colnames(df_bio)[3])
    }

    # Now, we generate the model:
    fit <- lmer(formula = paste(fixedFormula, "+", randomFormula), data = df_bio, REML = F)

    # And we extract the summary of the model:
    sum_table <- data.frame(coef(summary(fit)))

    # We change the name of the p-value column:
    colnames(sum_table)[5] <- "p-value"

    # We compute the R2 coefficient of the model and its components using rsq():
    R2 <- rsq(fit)

    # And then we format our summary table with the model information:
    sum_table$biomarker <- names(df)[i]
    sum_table$coefficient <- rownames(sum_table)
    if(!is.null(covariates)) {
      sum_table$covariate <- covariates
    }
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

  # Finally, we return the dataframe with all the information:
  return(lme_df)

}
