#' covBoxplot
#'
#' This function generates several boxplots that analyze the distribution of
#' numerical variables of a given dataframe with respect to the categorical
#' variables. It also performs parametric (Student's t or ANOVA depending
#' on the levels of the categorical variable) or non-parametric tests (Wilcoxon
#' or Kruskal-Wallis) to assess if the differences between the distributions
#' are significant. It shows in the boxplot's subtitle the name of the test
#' that has been used, the p-value and the sample size for each level of the
#' categorical variable. The function takes four arguments: the dataframe ("df"),
#' the range of columns with metadata information ("metaCols"), the option to
#' perform parametric or non-parametric tests ("parametric") and an
#' optional argument called "exportPDF" (NULL by default) that takes the
#' name of the PDF file that will store the plots. Finally, the function returns
#' a table with a summary of the distribution of numeric variables depending
#' on the levels of the categorical variables.
#'
#' @param df A given dataframe.
#' @param metaCols Range of columns with metadata information (e.g., 2:15).
#' @param parametric Option to perform parametric tests (Student's t, ANOVA)
#' or non-parametric (Wilcoxon or Kruskal-Wallis). TRUE by default.
#' @param exportPDF Export plots to a PDF (TRUE by default; if FALSE, plots will
#' be printed).
#' @param resultsFolder Name of the directory in which the results will be stored
#' ("results" by default).
#' @return A table with a summary of the distribution of numeric variables
#' depending on categorical variables
#' @importFrom tidyr pivot_wider
#' @import ggplot2 dplyr
#' @rawNamespace import(kableExtra, except = group_rows)
#' @export

covBoxplot <- function(df, metaCols, parametric = TRUE, exportPDF = TRUE, resultsFolder = "results") {
  # First, we check the existency of the resultsFolder to store the results:
  if (!dir.exists(resultsFolder)) {
    dir.create(resultsFolder)
  }

  # Then, we look at the exportPDF argument to see if we have to create a PDF
  # file to store all the plots:
  if(isTRUE(exportPDF)) {
    pdf(paste0(resultsFolder, "/covariates_boxplot.pdf"), height = 5, width = 10)
    cat("\n")
    cat(paste0("Plots were exported to the covariates_boxplot.pdf file (", resultsFolder, " folder)."))
    cat("\n\n")
  }

  # Then, we get the numeric variables and the categorical variables:
  cont_cov_names <- df[, metaCols] |>
    select(where(is.numeric)) |>
    colnames()

  # And then the categorical variables:
  cat_cov_names <-
    df[, metaCols] |>
    select(where(is.factor)) |>
    colnames()

  # Then we create a dataframe to store the results:
  statdf <- data.frame()

  # And also a common theme for the plots:
  bptheme <- theme(panel.background = element_rect(fill = "white"),
                   panel.grid.minor = element_blank(),
                   panel.grid.major = element_blank(),
                   axis.line = element_line(colour = "black", linewidth = 1),
                   plot.title = element_text(face = "bold"),
                   axis.text = element_text(size = 10),
                   axis.title = element_text(size = 12))

  # We set a color palette for the plots:
  plot_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")

  # Now, we iterate through every numeric variable:
  for (var1 in cont_cov_names) {
    for (var2 in cat_cov_names){
      df2 <- df |> filter(!is.na(df[var2]))
      n <- nrow(df2)
      # If there are two levels, we will do a Wilcoxon Test or a Student's t test
      # depending on the "parametric" argument:
      if(n_distinct(df2[var2]) == 2) {
        if(isTRUE(parametric)) {
          pval <- round(t.test(df2[,var1]~as.factor(df2[,var2]))$p.value,3)
          test_used <- "Student's T"
        } else {
          pval <- round(wilcox.test(df2[,var1]~as.factor(df2[,var2]))$p.value,3)
          test_used <- "Wilcoxon"
        }

        # If there are more than two levels, we will do an ANOVA or a Kruskal-Wallis
        # depending on the value of the "parametric" argument:
      } else if(n_distinct(df2[var2]) > 2) {
        if(isTRUE(parametric)) {
          anova_t <- aov(df2[,var1]~ as.factor(df2[, var2]))
          pval <- round(summary(anova_t)[[1]][["Pr(>F)"]][1], 3)
          test_used <- "ANOVA"
        } else {
          pval <- round(kruskal.test(df2[,var1]~ as.factor(df2[, var2]))$p.value, 3)
          test_used <- "Kruskal-Wallis"
        }

        # And if there is only one level, we will print a message and continue:
      } else {
        cat(paste("The variable", var2, "only has one unique level."))
        next
      }

      # We calculate the sample size for each level of the variable:
      counts <- df2 |>
        group_by(df2[, var2]) |>
        summarise(n = n())

      # And we create a vector to store the sample size for each level:
      n_labels <- paste0(paste0("N", 1:length(counts$n)), " = ", counts$n)

      # Now, we plot the results:
      p <- ggplot(data = df2, aes(df2[, var2], df2[, var1], color = df2[, var2]), na.rm = TRUE) +
        geom_boxplot(outlier.shape = NA, linewidth = 1) +
        geom_jitter(size = 1.5, alpha = 0.3, width = 0.3) +
        xlab(var2) +
        ylab(var1) +
        ggtitle(label = paste(var1, "distribution depending on", var2),
                subtitle = paste0(test_used, " test P value = ", pval, " (", paste(n_labels, collapse = ", "), ")")) +
        bptheme +
        scale_color_discrete(name = var2, type = plot_colors)



      # And we print the plot in the PDF file:
      print(p)

      # Lastly, we store variable, test, and P value data in the dataframe:
      statdf <- rbind(statdf, data.frame(var1 = var1, var2 = var2, test = test_used, p_value = pval))
    }
  }

  # We close the PDF file if we are exporting the plots:
  if(isTRUE(exportPDF)) {
    dev.off()
  }

  # Finally, we will use pivot_wider to obtain a more concise dataframe:
  results <- statdf |> pivot_wider(names_from = var1, values_from = p_value)

  # Then we will generate a Kable table to summarize our findings:
  tab <- results |> mutate(across(where(is.numeric),
                                    function(x){cell_spec(x,
                                                          color = ifelse(x < 0.05, "red", x),
                                                          bold = ifelse(x < 0.05, T, F))})) |>
    kable(escape = F, booktabs = T, caption = "Summary of the distribution of numeric variables depending on categorical variables",
          col.names = c("Categorical variable", "Statistical test", colnames(results[3:ncol(results)]))) |>
    kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = TRUE, position = "center") |>
    add_header_above(header = c(" " = 2, "P values" = length(cont_cov_names)))

  # And we return the table:
  return(tab)

}
