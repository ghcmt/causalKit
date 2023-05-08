#' covBarplot
#'
#' This function generates several bar plots that analyze the distribution of
#' categorical variables of a given dataframe with respect to the other
#' categorical variables. In addition, it performs a Chi-squared test to
#' compare if the differences between the two categorical variables tested
#' in each plot are significant or not. It provides the p-value of the test
#' as the subtitle, as well as the sample size for each group. It takes only
#' three arguments: the dataframe that we want to analyze ("df"), the range of
#' columns in which we have the covariates information ("metaCols") and an
#' optional argument called "exportPDF" (NULL by default) that takes the
#' name of the PDF file that will store the plots. Finally, the function returns
#' a table with a summary of the distribution of categorical variables depending
#' on the other categorical variables.
#'
#' @param df A provided dataframe.
#' @param metaCols Range of columns with metadata information (e.g., 2:15).
#' @param exportPDF Export plots to a PDF (TRUE by default; if FALSE, plots will
#' be printed).
#' @param resultsFolder Name of the directory in which the results will be stored
#' ("results" by default).
#' @return A table with a summary of the distribution of the categorical
#' variables of the dataframe depending on the other categorical variables.
#' @import ggplot2 dplyr
#' @rawNamespace import(kableExtra, except = group_rows)
#' @export

covBarplot <- function(df, metaCols, exportPDF = TRUE, resultsFolder = "results"){
  # First, we check the existency of the resultsFolder to store the results:
  if (!dir.exists(resultsFolder)) {
    dir.create(resultsFolder)
  }

  # Then, we look at the exportPDF argument to see if we have to create a PDF
  # file to store all the plots:
  if(isTRUE(exportPDF)) {
    pdf(paste0(resultsFolder, "/covariates_barplot.pdf"), height = 5, width = 10)
    cat("\n")
    cat(paste0("Plots were exported to the covariates_barplot.pdf file (", resultsFolder, " folder)."))
    cat("\n\n")
  }

  # We get only the categorical variables:
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

  # Palette for the plot:
  plot_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")


  # Now we iterate through each categorical variable:
  for (i in 1:length(cat_cov_names) ){
    # To compare with other variables while preventing the repetition of variables:
    if (i < length(cat_cov_names)){
      for (j in (i + 1):length(cat_cov_names)){
        # We remove NA values from the comparison:
        df2 <- df |> filter(!is.na(df[cat_cov_names[j]]) & !is.na(df[cat_cov_names[i]]))

        # We calculate the p-value with a chi-square test:
        pval <- round(chisq.test(df2[, cat_cov_names[i]], df2[, cat_cov_names[j]])$p.value, 3)

        # We calculate the sample size for each level of the variable:
        counts <- df2 |>
          group_by(df2[, cat_cov_names[i]]) |>
          summarise(n = n())

        # And we create a vector to store the sample size for each level:
        n_labels <- paste0(paste0("N", 1:length(counts$n)), " = ", counts$n)

        # We generate the plots:
        p <- ggplot(data = df2, aes(x = as.factor(df2[, cat_cov_names[i]]), fill = as.factor(df2[, cat_cov_names[j]]))) +
          geom_bar(position="dodge", color = "black", linewidth = 1) +
          ggtitle(label = paste(cat_cov_names[i], "distribution depending on", cat_cov_names[j]),
                  subtitle = paste0("Chi-squared test P value = ", pval, " (", paste(n_labels, collapse = ", "), ")")) +
          geom_text(aes(label= after_stat(count)), stat = 'count',
                    position = position_dodge(width = 0.9), vjust = -0.5, hjust = 0.5,
                    size = 4, family = "Helvetica", fontface = "bold") +
          bptheme + xlab(cat_cov_names[i]) +
          scale_fill_discrete(name = cat_cov_names[j], type = plot_colors)

        # And we print the plot:
        print(p)

        # Lastly, we store variables and P value data in the dataframe:
        statdf <- rbind(statdf, data.frame(var1 = cat_cov_names[i], var2 = cat_cov_names[j], p_value = pval))

      }
    }
  }

  # We close the PDF file if we are exporting the plots:
  if(isTRUE(exportPDF)) {
    dev.off()
  }

  # We order the dataframe according to ascendant p-value:
  statdf <- statdf |> arrange(p_value)

  # And we generate a table with the results:
  tab <- statdf |> mutate(across(where(is.numeric),
                                   function(x){cell_spec(x,
                                                         color = ifelse(x < 0.05, "red", x),
                                                         bold = ifelse(x < 0.05, T, F))})) |>
    kable(escape = F, booktabs = T, caption = "Summary of the distribution of categorical variables depending on other categorical variables") |>
    kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = TRUE, position = "center")

  return(tab)

}

