#' summarySigBio
#'
#' This function creates a table with the parameters of the model with and without
#' each covariate for selected biomarkers. The function requires the following
#' input parameters: "bio_names" (a vector with the names of the biomarkers to test),
#' "dfInt" (a dataframe that contains the values of the interaction linear models
#' (without covariates)), "dfMed" (a dataframe that contains the values of the
#' mediation linear models (with each covariate)), "ind_cov" (individual covariates
#' for which the interaction with the biomarkers is tested, and which will be
#' included as fixed effects in all models), "covariates" (a vector with the covariates
#' that are included in the models), "excel" (the name of the Excel file to store the
#' results), and tableCaption (the name of the table that returns the function).

#' The function iterates through the biomarkers and creates a subset of the interaction
#' dataframe with the target biomarkers. It then creates another dataframe from the
#' mediation dataframe with the interaction, and another from the mediation dataframe
#' with the information about covariates alone. Both dataframes are merged by the
#' target variables, and the level information is obtained for categorical variables.
#' The function then appends the information to the dataframe.

#' Next, the function calculates the change in percentage from the interaction estimate
#' without the covariate, and constructs the table, highlighting the rows with
#' an absolute change (comparing the estimate from the mediation dataframe
#' to the original estimate of the interaction dataframe without the covariate)
#' higher than 10% in red. The rows are grouped depending on the biomarker.
#' In addition, the function writes the results to an Excel file. Finally,
#' it returns a Kable table with the relevant information about the changes
#' in estimates.
#'
#' @param bio_names Vector with the names of the biomarkers to test.
#' @param dfInt Dataframe that contains the values of the interaction linear models
#' (without covariates).
#' @param dfMed Dataframe that contains the values of the mediation linear models
#' (with each covariate).
#' @param ind_cov Individual covariates for which we want to test the interaction
#' with the biomarkers. They will be included as fixed effects in all models.
#' @param covariates Vector with the covariates that we included in the models.
#' @param excel Name of the Excel file to store the results (e.g.,
#' "biomarkers_mediation_results.xlsx")
#' @param tableCaption Name of the table that returns the function (e.g., "Change
#' in Estimate value (%) after introducing covariates in the model for all
#' biomarkers").
#' @param resultsFolder Name of the directory in which the results will be stored
#' ("results" by default).
#' @return A Kable table with the parameters of the model with and without
#' each covariates for the selected biomarkers.
#' @import dplyr stringr openxlsx
#' @rawNamespace import(kableExtra, except = group_rows)
#' @export

# This function is designed to work inside LME_mediation:
summarySigBio <- function(bio_names, dfInt, dfMed, ind_cov, covariates, excel,
                          tableCaption, resultsFolder = "results") {
  # We create a dataframe to store the results:
  results <- data.frame()

  # We iterate through the biomarkers:
  for(bio in bio_names) {
    # We create a subset of the interaction dataframe with the target biomarkers:
    subset_1 <- dfInt |> filter(biomarker == bio)

    # Then, we create another from the mediation dataframe with the interaction:
    subset_2 <- dfMed |>
      filter(biomarker == bio & str_count(coefficient, ":") == length(ind_cov)) |>
      arrange(covariate)

    # And another dataframe from the mediation dataframe with the information
    # about covariates alone:
    subset_var <- dfMed |>
      filter(biomarker == bio & str_detect(coefficient, str_c(covariates, collapse = "|"))) |>
      arrange(covariate) |>
      rename(`p-value_var` = `p-value`, Estimate_var = Estimate, coef = coefficient) |>
      select(biomarker, covariate, Estimate_var, `p-value_var`, coef)

    # We merge both dataframes by the target variables. We set the "multiple"
    # argument as "all" as it is possible that we have several rows with the
    # same covariate value, as some covariates can have more than two levels.
    merge_df <- left_join(subset_2, subset_var, by = c("biomarker", "covariate"), multiple = "all")

    # We get the level information for categorical variables, as we want to know
    # from which level is the estimate and p-value information. If the variable
    # is numerical, the covariate value will be the same as the coef variable, so:
    merge_df$cov_level <- ifelse(merge_df$covariate != merge_df$coef, merge_df$coef, NA)

    # Next, we append the information to our dataframe. In the "covariate" column,
    # we will add between parentheses the information of the level if the variable
    # is categorical:
    results <- bind_rows(results, data.frame(biomarker = bio,
                                             estimate_int = subset_1$Estimate,
                                             p_val_int = subset_1$`p-value`,
                                             covariate = ifelse(is.na(merge_df$cov_level), merge_df$covariate,
                                                                paste0(merge_df$covariate, " (", merge_df$cov_level, ")")),
                                             estimate_cov = round(merge_df$Estimate_var, 3),
                                             p_val_cov = round(merge_df$`p-value_var`, 3),
                                             estimate_Int_Cov = round(merge_df$Estimate, 3),
                                             p_val_Int_Cov = round(merge_df$`p-value`, 3)))
  }

  # Now, we calculate the change in percentage from the interaction estimate
  # without the covariate:
  results$Change_pct <- round((results$estimate_Int_Cov - results$estimate_int)/(results$estimate_int)*100, 3)

  write.xlsx(results, file = paste0(resultsFolder, "/", excel))
  cat(paste0("Results were exported to the ", excel, " file (", resultsFolder, " folder)."))
  cat("\n\n")

  # Next, we begin to construct the table. We get the estimate and p-values information
  # about the biomarkers from the interaction dataframe (our starting point, as
  # in this dataframe we did not include covariates):
  ests <- results |> distinct(estimate_int) |> pull()
  pvals <- results |> distinct(p_val_int) |> pull()

  # We also are going to show in red those rows with an absolute change higher
  # than 10%:
  color_raw <- which(abs(results$Change_pct) > 10)

  # And we are going to group rows depending on the biomarker:
  packIndex <- rep(paste0(bio_names, " (Est. Int. = ", round(ests, 3), ", P-value Int. = ", round(pvals, 3), ")"), each = n_distinct(results$covariate))

  # We name the columns of the table:
  coltable <- c("Covariate", "Estimate Covariate", "P-value covariate",
                "Estimate Interaction", "P-value interaction", "Change (%)")

  # And we generate the table. We are going to show the information about the
  # first columns of the dataframe in the header of each group of rows:
  tab <- kable(results[, 4:ncol(results)], row.names = F,
               caption = tableCaption, col.names = coltable) |>
    kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                  full_width = TRUE, position = "center") |>
    row_spec(color_raw, bold = TRUE, color = "red") |>
    pack_rows(index = table(packIndex), background = "#fde9b7") |>
    add_header_above(header = c(" " = 1, "Covariate values" = 2, "Interaction + covariate values" = 2, " " = 1)) |>
    add_footnote("The level of the variable is specified between parentheses if the variable is categorical.", notation = "symbol")

  # Finally, we return the table:
  return(tab)
}
