#' covAnalysis
#'
#' This function analyzes the influence of the covariates of a given dataframe
#' in the levels of the selected biomarkers using linear mixed effects models.
#' The main purpose of this function is to detect the existence of potential
#' confounding covariates. It takes several arguments as input: a given
#' dataframe ("df"), the range of columns with the covariates (metaCols),
#' the name of the patient or ID column (patientCol), the range of columns
#' that contain the biomarker data (bioCols), the individual covariates which
#' will be included as an interaction the model (ind_cov), the method used to
#' evaluate the models (REML) and an optional argument (TRUE by default) that
#' will print the formula of each model.
#'
#' The function iterates through each metadata variable (potential confounding
#' variables), generates the formulas for the linear mixed-effects models, and
#' calls the lmeBioCov function that will generate a model for all biomarkers.
#' The results from each iteration are added to the results dataframe. Finally,
#' the results are exported to an Excel file and the function returns the
#' a dataframe with all the models generated.
#'
#' @param df A given dataframe.
#' @param metaCols Range of columns with metadata information (e.g., 2:15)
#' @param patientCol Name of the ID or patient column in the top-table.
#' @param bioCols Range of columns with biomarker information (e.g., 16:300)
#' @param ind_cov Individual covariates for which we want to test the interaction
#' with the biomarkers. They will be included as fixed effects in all models.
#' @param REML Algorithm used for estimating the model parameters. It is recommended
#' to set it as FALSE (it will use ML) if we are comparing models.
#' @param printFormula Print the formula of the different models (TRUE by default).
#' @param resultsFolder Name of the directory in which the results will be stored
#' ("results" by default).
#' @return A dataframe with the results of the mediation analysis
#' @import dplyr openxlsx lmerTest
#' @export


covAnalysis <- function(df, metaCols, patientCol, bioCols, ind_cov, REML, printFormula = TRUE,
                        resultsFolder = "results") {
  # First, we get the names of the proteins variables from the bioCols:
  bio_names <- colnames(df[bioCols])

  # And the names of the metadata variables:
  meta_vars <- colnames(df[metaCols])

  # We remove the individual covariates given by the user from those vars
  # (in case that they are there):
  meta_vars <- meta_vars[!(meta_vars %in% ind_cov)]

  # We also check the existency of the resultsFolder to store the results:
  if (!dir.exists(resultsFolder)) {
    dir.create(resultsFolder)
  }

  # Then, we create a dataframe to store the information:
  results_df <- data.frame()

  # We iterate through every variable:
  for(var in meta_vars) {
    # We get the index of the variable:
    i <- which(names(df) == var)

    # Now we generate the formulas of the model:
    formula1 <- paste("biomarker ~", paste0(ind_cov, collapse = "*"), "+", var)
    formula2 <- paste0("(1 | ", patientCol, ")")

    # We show the formula of the model if desired:
    if(isTRUE(printFormula)) {
      message(paste0("\n ", formula1, " + ",  formula2, " \n "))
    }


    # And we call the LME_results_P1 function. It is important to omit NA
    # values as they can cause an error:
    lme_df <- lmeBioCov(df = df,
                        bio_names = bio_names,
                        ind_cov = ind_cov,
                        covariate = var,
                        patientCol = patientCol,
                        fixedFormula = formula1,
                        randomFormula = formula2,
                        REML = REML)

    # And we add the newly generated rows to the general dataframe:
    results_df <- rbind(lme_df, results_df)

  }

  results_df <- results_df |> relocate(coefficient, .before = 1)

  # Then, we export the results to an Excel file:
  cat("\n")
  write.xlsx(results_df, file = paste0(resultsFolder, "/confounding_covariates_results.xlsx"))
  cat(paste0("Full confounding analysis results are stored in the confounding_covariates_results.xlsx file (", resultsFolder, " folder)."))
  cat("\n\n")

  # We generate a table with the first rows of the analysis:
  modelTable <- kable(results_df[1:10, ], row.names = F,
                      caption = "First 10 rows of the Confounding Analysis") |>
    kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                  full_width = TRUE, position = "center")

  # And then we print it:
  print(modelTable)

  # Finally, we return the dataframe with all the data:
  return(results_df)
}
