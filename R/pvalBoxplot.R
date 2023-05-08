#' pvalBoxplot
#'
#' This function generates a boxplot and a table summarizing the number of
#' significant biomarkers by covariate after running a linear model. Besides the
#' "df", it takes as input a vector of "covariates" (the vector with the covariates
#' that we will put into the model) "pvalCol" (indicating the name of the p-value
#' column), and a "covCol" (which indicates the name of the Column containing the
#' covariate information). The function also has an optional argument, "exportPDF" (which specifies the file name to which
#' the plot should be exported as a PDF file), and "sigs" (which is a vector of
#' significance levels to be used in the output table).

#' @param df A dataframe with the results of a linear model.
#' @param covariates Vector with the covariates that we want to include in the model.
#' @param pvalCol Name of the p-value column in the dataframe.
#' @param covCol Name of the covariates column in the dataframe.
#' @param exportPDF Name of the PDF file in which the boxplot is exported to
#' (NULL by default; if NULL, it will print the plot).
#' @param sigs Vector with the desired adjusted p-value significance levels
#' (e.g., c(0.01, 0.05)).
#' @param resultsFolder Name of the directory in which the results will be stored
#' ("results" by default).
#' @return A table with the number of biomarkers below the desired significance
#' thresholds (adjusted p-value) for each covariate.
#' @importFrom tidyr pivot_wider
#' @import ggplot2 dplyr
#' @rawNamespace import(kableExtra, except = group_rows)
#' @export


pvalBoxplot <- function(df, covariates, pvalCol, covCol, exportPDF = TRUE, sigs, resultsFolder = "results") {
  # First, we create the dataframe to store the data:
  adj_df <- data.frame()

  # We also check the existency of the resultsFolder to store the results:
  if (!dir.exists(resultsFolder)) {
    dir.create(resultsFolder)
  }

  # Then, we generate the subset with the desired covariates:
  for (i in covariates) {
    subset_df <- df |> filter(str_detect(.data[[covCol]], paste0(i, ".*?")))
    subset_df$adj_p_val <- p.adjust(subset_df[[pvalCol]], method = "BH")
    adj_df <- bind_rows(adj_df, subset_df)
  }

  if(isTRUE(exportPDF)) {
    pdf(paste0(resultsFolder, "/significance_boxplot.pdf"), height = 5, width = 10)
    cat("\n")
    cat(paste0("Plots were exported to the significance_boxplot.pdf file (", resultsFolder, " folder)."))
    cat("\n\n")
  }

  colorPalette <- c("#5C4B51", "#8CBEB2", "#8CBEB2", "#F3B562", "#F06060",
               "#a8949b", "#bedad3", "#f8f4db", "#f6c98c", "#f48a8a",
               "#60f0f0", "#4b5c56", "#be8c98", "#bfc6f2", "#62a0f3",
               "#60f0a8", "#5a5c4b", "#b1be8c", "#bff2d2", "#62f36d",
               "#f0f060", "#A13323", "#2391a1", "#a346fb", "#E2A72E")

  # We generate the plot:
  p1 <- ggplot(adj_df, aes(x = .data[[covCol]], y = -log(.data[[pvalCol]], 10))) +
    geom_boxplot(aes(fill = .data[[covCol]]), linewidth = 0.5) +
    theme_minimal() +
    scale_fill_manual(values = colorPalette) +
    geom_hline(yintercept = c(4, 5),lty = 2, col = c("blue","red")) +
    theme(axis.text = element_text(size = 15),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.title = element_text(size = 15, face = "bold"),
          legend.position = "none",
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "white"),
          axis.line = element_line(colour = "black", linewidth = 1)) +
    xlab("covariates") + ylab("significance (-log10(p.value)") +
    ggtitle("P-value distribution for each coefficient",
            subtitle = "P-value below 0.0001 (blue line) or 0.00001 (red line)")

  print(p1)

  # We close the PDF file if we are exporting the plots:
  if(isTRUE(exportPDF)) {
    dev.off()
  }

  # Finally, we create a dataframe to store the significance information:
  sig_df <- data.frame()

  # And we iterate through the significance levels:
  for (sig in sigs) {
    int_df <- adj_df |>
      group_by(coefficient) |>
      summarise(significance = sig,
                number = sum(adj_p_val < sig),
                n = n())
    sig_df <- bind_rows(sig_df, int_df)
  }

  # We change the format of the dataframe to have one column for each significance
  # level introduced. We move the sample column to the end, as we don't want
  # to have an additional column in the table with this information:
  sig_df <- pivot_wider(sig_df, names_from = significance, values_from = number) |> relocate(n, .after = last_col())

  # And we create the table. We show all the columns except for the last one, as we are
  # going to show the information about the sample size in the title of the table.
  sigTable <- kable(sig_df[, -ncol(sig_df)], row.names = F,
                    caption = paste0("Number of changed biomarkers by significance level (adjusted p-value) (N = ",
                                     max(sig_df$n), ")"),
                    col.names = c("Coefficient", sigs)) |>
    kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = TRUE, position = "center") |>
    add_header_above(header = c(" " = 1, "Level of significance (adj. p-val)" = length(sigs)))

   # We print the table:
  print(sigTable)

  # Finally, we return the dataframe:
  return(sig_df)

}
