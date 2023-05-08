#' modelCompLRT
#'
#' This function takes in several arguments: a dataframe "df", a vector of
#' covariates ("covariates"), the name of the ID or patient column ("patientCol"),
#' a range of columns with biomarker information ("bioCols"), individual covariates
#' for which we want to test the interaction with the biomarkers ("ind_cov"), the
#' algorithm used for estimating the model parameters ("method"), and a flag to
#' indicate whether to print the formula of the different models ("printFormula").
#' The purpose of this function is to find if the interaction of some covariates
#' with the individual covariates has an influence in the dependent variable,
#' which in this case is a list of biomarkers.

#' For each covariate, the function generates two models: one
#' including the interaction between the covariate and the individual covariates
#' and another without the interaction. The function then calls the lmeModelComp
#' function to compare the models, and stores the results in a dataframe.
#' The function then summarizes the number of biomarkers with a p-value below
#' 0.05 for each covariate, and returns a table with this summary.
#'
#' @param df A given dataframe
#' @param covariates Vector with the covariates that we want to include in the model.
#' @param patientCol Name of the ID or patient column in the top-table.
#' @param bioCols Range of columns with biomarker information (e.g., 16:300)
#' @param ind_cov Individual covariates for which we want to test the interaction
#' with the biomarkers. They will be included as fixed effects in all models.
#' @param REML Algorithm used for estimating the model parameters. It is recommended
#' to set it as FALSE (it will use ML) if we are comparing models.
#' @param printFormula Print the formula of the different models in the console
#' (TRUE by default).
#' @param resultsFolder Name of the directory in which the results will be stored
#' ("results" by default).
#' @return A table with the summary of the model comparison.
#' @importFrom dplyr group_by summarise
#' @import openxlsx lmerTest
#' @rawNamespace import(kableExtra, except = group_rows)
#' @export


modelCompLRT <- function(df, covariates, patientCol, bioCols, ind_cov, REML,
                         printFormula = TRUE, resultsFolder = "results") {
  # First, we get the names of the proteins variables from the bioCols:
  bio_names <- colnames(df[bioCols])

  # Then, we create a dataframe to store the information:
  results_df <- data.frame()

  # We also check the existency of the resultsFolder to store the results:
  if (!dir.exists(resultsFolder)) {
    dir.create(resultsFolder)
  }

  # And we initiate a counter to change the formula depending on the coefficient:
  i <- 1

  # We iterate through every variable:
  for(cov in covariates) {
    # We store the rest of coefficients in a new variable:
    other_cov <- covariates[-i]

    # Now we generate the two formulas:
    formula1 <- paste0("biomarker ~ ", paste0(ind_cov, collapse = "*"), "*", cov, "+", paste0(other_cov, collapse = " + "))
    formula2 <- paste0("biomarker ~ ", paste0(ind_cov, collapse = "*"), "+", cov, "+", paste0(other_cov, collapse = " + "))

    # And we use the same random formula for all cases:
    rndForm <- paste0("(1 | ", patientCol, ")")

    # We show the formula of the model if desired:
    if(isTRUE(printFormula)) {
      message(paste0("\n ", formula1, " + ",  rndForm, " \n "))
      message(paste0("\n ", formula2, " + ",  rndForm, " \n "))
    }


    # And we call the lme_comp function. It is important to omit NA
    # values as they can cause an error:
    comp <- lmeBioLRT(df =  df,
                     bio_names = bio_names,
                     ind_cov = c(covariates, ind_cov),
                     patientCol = patientCol,
                     fixedFormula1 = formula1,
                     fixedFormula2 = formula2,
                     covariates = cov,
                     randomFormula = rndForm,
                     REML = REML)

    # And we add the newly generated rows to the general dataframe:
    results_df <- rbind(comp, results_df)

    # And we add to the counter:
    i <- i + 1

  }

  # We export to an Excel file the results:
  cat("\n")
  write.xlsx(results_df, file = paste0(resultsFolder, "/model_comparison_LRT.xlsx"))
  cat(paste0("Model comparison results are stored in the model_comparison_LRT.xlsx file (", resultsFolder, " folder)."))
  cat("\n\n")

  # With the final dataframe, we are going to create a table that summarizes
  # the number of biomarkers with a p-val below 0.05 for each coefficient:
  comp_summary <- results_df |>
    group_by(coefficient) |>
    summarise(number = sum(p_val < 0.05),
              n = n())

  # And we format the summary:
  compTable <- kable(comp_summary[, -ncol(comp_summary)], row.names = F,
                     caption = paste0("Number of biomarkers by coefficient with a p-value below 0.05 after comparing the model with or without the interaction (N = ", max(comp_summary$n), ")"), col.names = c("Coefficient", "Number of biomarkers")) |>
    kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = TRUE, position = "center")

  # And then we print it:
  print(compTable)

  # Finally, we return the dataframe with the model information:
  return(results_df)

}
