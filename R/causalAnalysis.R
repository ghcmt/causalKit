#' causalAnalysis
#'
#' This function runs all the functions in the package.
#'
#' @param df A given dataframe
#' @param bioCols Range of columns with biomarker information (e.g., 16:300)
#' @param ind_cov Individual covariates for which we want to test the interaction
#' with the biomarkers. They will be included as fixed effects in all models.
#' @param metaCols Range of columns with metadata information (e.g., 2:15)
#' @param patientCol Name of the ID or patient column in the top-table.
#' @param depVar Dependent variable, as clinical outcome or score (NULL by default)
#' @param colorVar A categorical variable from the dataframe to color the points.
#' (e.g., "Group"). NULL by default.
#' @param REML Algorithm used for estimating the model parameters. It is recommended
#' to set it as FALSE (it will use ML) if we are comparing models.
#' @param parametric Option to perform parametric tests (Student's t, ANOVA)
#' or non-parametric (Wilcoxon or Kruskal-Wallis) in boxplots. TRUE by default.
#' @param exportPDF Export plots to a PDF (TRUE by default; if FALSE, plots will
#' be printed).
#' @param printMedTable Print the table with the information about the biomarkers
#' that lost significance after adding a covariate (mediation analysis) (TRUE
#' by default).
#' @param resultsFolder Name of the directory in which the results will be stored
#' ("results" by default).
#' @return A dataframe with the results of the different linear models.
#' @import dplyr progress
#' @rawNamespace import(kableExtra, except = group_rows)
#' @export
#'

causalAnalysis <- function(df, bioCols, ind_cov, metaCols, patientCol, depVar = NULL,
                           colorVar = NULL, REML = FALSE, parametric = TRUE, exportPDF = TRUE,
                           printMedTable = TRUE, resultsFolder = "results") {

  # We start by creating a results folder if there isn't any:
  if (!dir.exists(resultsFolder)) {
    dir.create(resultsFolder)
  }

  cat("\n")
  cat("# Causal Analysis")
  cat("\n")

  # First, we show the correlation between continous variables:
  cat("## Correlation between numeric variables")
  cat("\n")
  p <- covPairs(df, metaCols, colorVar)
  print(p, out.width = 15, out.height = 15)
  cat("\n\n")

    # We ask the user if the covBoxplot function should be ran:
  message("The covBoxplot function generates boxplots with the distribution of continous
variables depending on categorical variables.")
  run_covBoxplot <- tolower(readline(prompt = "Do you want to run the covBoxplot function? (y/n): "))



  if (run_covBoxplot != "n") {
    # Now, the boxplots with the distribution of numeric variables depending
    # on categorical variables with stats on them:
    cat("## Distribution of continuous variables depending on categorical variables")
    cat("\n")
    covBoxplot(df, metaCols, parametric, exportPDF, resultsFolder)
    cat("\n")
  }

  # We ask the user if the covBarplot function should be ran:
  message("The covBarplot function generates bar plots with the distribution of categorical
variables depending on the other categorical variables.")
  run_covBarplot <- tolower(readline(prompt = "Do you want to run the covBarplot function?(y/n): "))

  if (run_covBarplot != "n") {
    # Next, barplots with the distribution of categorical variables depending
    # on the other categorical variables (and a Chi-squared test):
    cat("## Distribution of categorical variables depending on other categorical variables")
    cat("\n")
    covBarplot(df, metaCols, exportPDF, resultsFolder)
    cat("\n")
  }

  # We ask the user if the Confounding Analysis function should be performed:
  message("The covAnalysis function analyzes the effect on biomarkers of adding each covariate to the
model with the pre-defined individual covariates. It also will generate a boxplot with the results.")
  run_covAnalysis <- tolower(readline(prompt = "Do you want to run the covAnalysis function? (y/n): "))

  if (run_covAnalysis != "n") {
    # We pivot to mixed models. First, we are going to perform a covariate analysis
    # to establish which covariates should go into the model:
    cat("## Confounding Covariates Analysis")
    cat("\n")
    cat("### Mixed Models")
    bio_names <- colnames(df[bioCols])

    # We store the result of the function in a dataframe, as this will be the
    # input for the pvalBoxplot() function:
    lme_df <- covAnalysis(df, metaCols, patientCol, bioCols, ind_cov, REML,
                          printFormula = TRUE, resultsFolder)

    # We extract the covariates to give as an input to the following function::
    meta_vars <- colnames(df[metaCols])
    meta_vars <- meta_vars[!(meta_vars %in% ind_cov)]

    cat("\n")
    cat("### Significance boxplot")
    cat("\n")
    sig_df <- pvalBoxplot(df = lme_df, covariates = meta_vars, pvalCol = "p-value",
                          covCol = "coefficient", exportPDF, sigs = c(0.01, 0.05), resultsFolder)

    # We get the covariates that have more than 10% of biomarkers with an
    # adjusted p-value lower than the last significance column:
    covs <- sig_df |> filter((sig_df[, ncol(sig_df) - 1])/n > 0.1)
    covs <- str_extract(covs$coefficient, paste(meta_vars, collapse = "|"))
    cat("\n")
    cat(paste("The covariates which have more than a 10% of the total sample size with a p-value below 0.05
            are recommended to be included in the model. These variables are:", paste(covs, collapse = ", "), "."))
    cat("\n")
  } else {
    covs <- colnames(df[metaCols])
  }

  # We ask the user if the Model Comparison should be performed:
  message("The modelCompLRT function performs Model Comparison Analysis. It will create two models by
adding each covariate to the model: one inside the interaction and one without
the interaction. Then, it will compare the two models using a Likelihood Ratio Test")
  run_modelComp <- tolower(readline(prompt = "Do you want to run the modelCompLRT function? (y/n): "))

  if (run_modelComp != "n") {
    # Now, we are going to compare models with and without the covariates as
    # an interaction term:
    cat("\n")
    cat("## Model Comparison")
    cat("\n")
    modelCompLRT(df, covariates = covs, patientCol, bioCols, ind_cov, REML,
                 printFormula = TRUE, resultsFolder)
    cat("\n\n")
  }



  # Finally, we will perform the mediation analysis if there is a dependent variable:
  if(!is.null(depVar)) {
    # We ask the user if the analysis should be ran:
    message("The mediationAnalysis function analyzes the effect of all biomarkers in the dependent variable
by adding them to the model inside the main interaction. Then, it will look at how many biomarkers
loss significance after adding a covariate to the model.")
    run_mediationAnalysis <- tolower(readline(prompt = "Do you want to run the mediationAnalysis function? (y/n): "))

    if (run_mediationAnalysis != "n") {
      cat("\n\n")
      cat("## Mediation Analysis")
      cat("\n")
      # First, we get the names of the proteins variables from the bioCols:
      bio_names <- colnames(df[bioCols])
      med_df <- mediationAnalysis(df, metaCols, patientCol, bio_names, ind_cov,
                                  depVar, REML, printFormula = TRUE,
                                  printMedTable, resultsFolder)
    }
  }

}
