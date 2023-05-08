#' lmeBioLRT
#'
#' This function performs a likelihood ratio test (LRT) to compare the fit of two
#' linear mixed-effect models for multiple biomarkers in a given dataframe. The
#' function uses the biomarkers as the dependent variable and it compares
#' the differences between a model that adds a covariate to a given interaction
#' and the model with the covariate added to the model but outside of the
#' main interaction. The function takes in a dataframe "df" and several other
#' parameters, including "bio_names" (a vector with the names of the biomarkers to
#' test), "covariates" (a vector with the covariates to include in the model),
#' "patientCol" (the name of the ID or patient column in the df), "ind_cov"
#' (the individual covariates for which the interaction with the biomarkers will
#' be tested and included as fixed effects in all models), "fixedFormula1" and
#' "fixedFormula2" (the formulas for the first and second fixed effects of the
#' model), "randomFormula" (the formula for the random effects of both
#' models), and "REML" (the algorithm used for estimating the model parameters).

#' The function first retrieves the indexes of the biomarkers, patient, and
#' individual covariates and creates a data frame to store the data. It loops
#' through all the biomarkers, generates two models for each one, performs the
#' LRT, and stores the relevant data in a data frame. The function exports the
#' results to an Excel file and returns the data frame with the comparisons.
#' The function uses the lme function from the nlme package to fit the
#' mixed-effects models.
#'
#' @param df A given dataframe
#' @param bio_names Vector with the names of the biomarkers to test.
#' @param covariates Vector with the covariates that we want to include in the model.
#' @param patientCol Name of the ID or patient column in the dataframe.
#' @param ind_cov Individual covariates for which we want to test the interaction
#' with the biomarkers. They will be included as fixed effects in all models.
#' @param fixedFormula1 Formula of the first Fixed Effects of the Model
#' (with the covariate added to the interaction).
#' @param fixedFormula2 Formula of the second Fixed Effects of the Model
#' (without the covariate added to the interaction).
#' @param randomFormula Formula of the Random Effects of both Models.
#' @param REML Algorithm used for estimating the model parameters. It is recommended
#' to set it as FALSE (it will use ML) if we are comparing models.
#' @return A dataframe with the summary of the model comparison.
#' @import openxlsx lmerTest rsq progress
#' @importFrom nlme fixed.effects
#' @export


lmeBioLRT <- function(df, bio_names, covariates, patientCol, ind_cov, fixedFormula1,
                      fixedFormula2, randomFormula, REML) {
  # First, we get the indexes of the biomarkers, patient and individual covariates:
  vars <- names(df)
  bio_id <- which(vars %in% bio_names)
  patient_id <- which(names(df) == patientCol)
  cov_id <- which(vars %in% ind_cov)

  # And we create a dataframe to store the data:
  comp_df <- data.frame()

  # We start a counter:
  j <- 0

  # We create a progress bar to show the time remaining for every model:
  pb <- progress_bar$new(total = length(bio_id), format = "[:bar] :percent :eta")

  # And we loop through all of the biomarkers:
  for (i in bio_id){
    # First, we add to the counter:
    j <- j + 1

    # We create a subset with the target columns: biomarker, patient and covariates:
    df_bio <- df[, c(patient_id, i, cov_id)]

    # We create a dataframe without NA to calculate the R2 coefficient:
    df_clean <- na.omit(df_bio)

    # We change the name of the biomarker column to "biomarker" to have the same
    # name as in the model formula:
    colnames(df_bio)[2] <- "biomarker"

    # Now, we generate the model 1:
    fit1 <- lmer(formula = paste(fixedFormula1, "+", randomFormula), data = df_bio, REML = F)

    # And model 2:
    fit2 <- lmer(formula = paste(fixedFormula2, "+", randomFormula), data = df_bio, REML = F)

    # We perform the calculations and we store the values in the dataframe. First,
    # we calculate the degrees of freedom:
    deg_freedom <- length(fixed.effects(fit1)) - length(fixed.effects(fit2))

    # And the p-value using the logLik values from each model and the degrees
    # of freedom that we have just calculated:
    p_val <- pchisq(summary(fit1)$logLik[1] - summary(fit2)$logLik[1], df = deg_freedom, lower.tail = F, log.p = F)

    # Then, we store the relevant data in a dataframe:
    sum_df <- data.frame(
      deg_freedom = deg_freedom,
      p_val = p_val,
      biomarker = names(df)[i],
      index = bio_id[j],
      coefficient = covariates
    )

    # We bind the newly generated rows to the existent dataframe:
    comp_df <- bind_rows(comp_df, sum_df)

    # And we add a tick to the progress bar:
    pb$tick()

  }

  # Finally, we return the dataframe with the comparisons:
  return(comp_df)

}
