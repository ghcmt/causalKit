#' mediationAnalysis
#'
#' This function mediationAnalysis performs a mediation analysis between the biomarkers
#' and the dependent variable while considering the effect of individual covariates
#' on this relationship. It takes as input a data frame ("df"), the names of the
#' columns corresponding to metadata ("metaCols"), the column name containing the
#' patients' IDs ("patientCol"), the names of the columns corresponding to
#' biomarkers ("bioCols"), the individual covariates ("ind_cov"), the dependent
#' variable ("depVar"), the method used to perform the analysis ("REML"),
#' and a logical argument to indicate whether the formula of the model should be
#' printed ("printFormula", TRUE by default).
#'
#' First, the function identifies the biomarkers with significant interactions with
#' the individual covariates, and then, it performs again the same analysis but
#' adding the other covariates that the user wants to test. The results are stored
#' in a data frame and written to an Excel file in the "results" folder. The
#' function returns a dataframe with the biomarkers which interaction with the individual
#' covariates is significant in a model with the dependent variable.
#'
#' @param df A given dataframe
#' @param metaCols Range of columns with metadata information (e.g., 2:15)
#' @param patientCol Name of the ID or patient column in the top-table.
#' @param bio_names Vector with the names of the biomarkers to test.
#' @param ind_cov Individual covariates for which we want to test the interaction
#' with the biomarkers. They will be included as fixed effects in all models.
#' @param depVar Dependent variable, as clinical outcome or score.
#' @param REML Algorithm used for estimating the model parameters. It is recommended
#' to set it as FALSE (it will use ML) if we are comparing models.
#' @param printFormula Print the formula of the different models in the console
#' (TRUE by default).
#' @param printMedTable Print the table with the information about the biomarkers
#' that lost significance after adding a covariate (mediation analysis) (TRUE
#' by default).
#' @param resultsFolder Name of the directory in which the results will be stored
#' ("results" by default).
#' @return A dataframe with the results of the mediation analysis
#' @import dplyr stringr openxlsx lmerTest rsq
#' @rawNamespace import(kableExtra, except = group_rows)
#' @export

mediationAnalysis <- function(df, metaCols, patientCol, bio_names, ind_cov, depVar,
                              REML, printFormula = TRUE, printMedTable = TRUE,
                              resultsFolder = "results") {
  # We get the names of the metadata variables:
  meta_vars <- colnames(df[metaCols])

  # We remove the individual covariates given by the user from those vars
  # (in case that they are there):
  meta_vars <- meta_vars[!(meta_vars %in% ind_cov)]

  # We also check the existency of the resultsFolder to store the results:
  if (!dir.exists(resultsFolder)) {
    dir.create(resultsFolder)
  }

  # Next, we get the information about the interaction of our individual covariates
  # with the biomarkers in the dependent variable:
  formula1 <- paste(depVar, "~", paste0(ind_cov, collapse = "*"), "* biomarker")
  formula2 <- paste0("(1 | ", patientCol, ")")
  int_df <- lmeBioInt(df =  df,
                      bio_names = bio_names,
                      ind_cov = ind_cov,
                      patientCol = patientCol,
                      depVar = depVar,
                      fixedFormula = formula1,
                      randomFormula = formula2,
                      REML = REML)


  # We filter the data as we want to only retain the triple interaction:
  int_df <- int_df |> filter(str_count(coefficient, ":") == length(ind_cov))

  # Then we arrange by p-value:
  int_df <- int_df |>
    arrange(`p-value`) |>
    relocate(coefficient, .before = 1)

  # We export the results to an Excel file:
  cat("\n")
  write.xlsx(int_df, file = paste0(resultsFolder, "/biomarkers_interaction_analysis_", depVar, ".xlsx"))
  cat(paste0("The results about the biomarkers that have a significant interaction with individual covariates
  are stored in the ", paste0("/biomarkers_interaction_analysis_", depVar, ".xlsx"), " file (", resultsFolder, " folder)."))
  cat("\n\n")

  # And we gather the information in a table, coloring depending on p-value:
  color_raw <- which(int_df$`p-value` < 0.05)
  color_raw <- color_raw[color_raw <= 10]
  t <- kable(int_df[1:10, ], row.names = F, caption = paste0("Biomarkers with significant interaction with individual covariates")) |>
    kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = TRUE, position = "center") |>
    row_spec(color_raw, bold = TRUE, color = "red") |>
    footnote(paste0("Ordered by ascendent p-value."))

  # We print the table:
  print(t)

  # Now, we are going to perform the mediation analysis. First, we will get
  # the biomarkers with significant interaction to iterate through them:
  sig_bio <- int_df |> filter(`p-value` < 0.05) |> distinct(biomarker) |> pull()

  # Then, we create a dataframe to store the information of the analysis:
  results_df <- data.frame()

  # We iterate through every variable:
  for(var in meta_vars) {
    # We get the index of the variable:
    i <- which(names(df) == var)

    # Now we generate the formulas of the model:
    formula1 <- paste(depVar, "~", paste0(ind_cov, collapse = "*"), "* biomarker +", var)
    formula2 <- paste0("(1 | ", patientCol, ")")

    # We show the formula of the model if desired:
    if(isTRUE(printFormula)) {
      message(paste0("\n ", formula1, " + ",  formula2, " \n "))
    }

    # And we call the LME_bio_int function. It is important to omit NA
    # values as they can cause an error. In this case, we specify the coef
    # argument as we want to include covariates in the model:
    lme_df <- lmeBioInt(df = df,
                        bio_names = sig_bio,
                        ind_cov = ind_cov,
                        patientCol = patientCol,
                        depVar = depVar,
                        covariates = var,
                        fixedFormula = formula1,
                        randomFormula = formula2,
                        REML = REML)

    # And we add the newly generated rows to the general dataframe:
    results_df <- rbind(lme_df, results_df)

  }

  # We filter the data as we want to retain the values of the interaction and the values
  # of the covariates in the model:
  results_df <- results_df |>
    filter(str_count(coefficient, ":") == length(ind_cov) |
             str_detect(coefficient, str_c(meta_vars, collapse = "|"))) |>
    relocate(coefficient, .before = 1)

  # Then we arrange by p-value:
  results_df <- results_df |>
    arrange(`p-value`)

  # We export the results to an Excel file:
  cat("\n")
  write.xlsx(results_df, file = paste0(resultsFolder, "/mediation_analysis_significant_biomarkers.xlsx"))
  cat(paste0("The results about the biomarkers that mantain a significant interaction after adding individual covariates to the model
  are stored in the mediation_analysis_significant_biomarkers.xlsx file (", resultsFolder, " folder)."))
  cat("\n\n")

  # We color the rows with significant values:
  color_raw <- which(results_df$`p-value` < 0.05)
  color_raw <- color_raw[color_raw <= 10]

  # And we print the table:
  tab <- kable(results_df[1:10, ], row.names = F, caption = paste0("Biomarkers with significant interaction after adding individual covariates to the model")) |>
    kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = TRUE, position = "center") |>
    row_spec(color_raw, bold = TRUE, color = "red") |>
    footnote(paste0("Ordered by ascendent p-value."))
  print(tab)

  # Then, we assess the number of significant biomarkers depending on covariate.
  # We want to know how many biomarkers are still significant after adding each
  # covariate to the model:
  med_summary <- results_df |>
    filter(str_count(coefficient, ":") == length(ind_cov)) |>
    group_by(covariate) |>
    summarise(n = n(),
              num = sum(`p-value` < 0.05))

  # And we generate the table to print it:
  tab2 <- kable(med_summary, row.names = F, caption = "Significative biomarkers after introducing each covariate in the model",
                col.names = c("Covariate", "Total number", "P-val < 0.05")) |>
    kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = TRUE, position = "center")
  print(tab2)

  # Next, we get the non-significant biomarkers (i.e., those who lost significance
  # after introducing the covariate)
  non_sig_bio <- results_df |>
    filter(str_count(coefficient, ":") == length(ind_cov)) |>
    filter(`p-value` > 0.05) |>
    distinct(biomarker) |>
    pull() |>
    sort()

  # Then, we create a table showing the biomarkers in which the interaction
  # is not significant after adding covariates:
  if(length(non_sig_bio) > 0) {
    non_sig_table <- summarySigBio(bio_names = non_sig_bio,
                                   excel = "estimate_change_nonsignificant_biomarkers.xlsx",
                                   tableCaption = "Change in Estimate value (%) after introducing covariates in the model for biomarkers in which interaction significance was lost",
                                   dfInt = int_df,
                                   dfMed = results_df,
                                   ind_cov = ind_cov,
                                   covariates = meta_vars,
                                   resultsFolder = resultsFolder)

    # We show the table if desired:
    if(isTRUE(printMedTable)) {
      print(non_sig_table)
    }

  } else {
    print("No biomarker has lost significance when adding a covariate to the model.")
  }

  # And we also generate an Excel file with the information for all biomarkers:
  sig_bio <- sig_bio |> sort()

  if(length(sig_bio) > 0) {
    summarySigBio(bio_names = sig_bio,
                  excel = "estimate_change_all_biomarkers.xlsx",
                  tableCaption = "Change in Estimate value (%) after introducing covariates in the model for all biomarkers",
                  dfInt = int_df,
                  dfMed = results_df,
                  ind_cov = ind_cov,
                  covariates = meta_vars,
                  resultsFolder = resultsFolder)
  } else {
    print("No biomarkers have significant interaction with individual covariates.")
  }

  # Finally, we return the dataframe with the mediation results:
  return(results_df)
}
