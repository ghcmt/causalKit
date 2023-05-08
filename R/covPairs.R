#' covPairs
#'
#' This function is an adaptation of the ggpairs function from the GGally package
#' that eases its use with any given dataframe. It needs a dataframe, the range
#' of the metadata columns and an optional grouping variable to color the points
#' and plots. From the range of the metadata columns, the function extracts
#' the numeric variables and generates the plot.
#'
#' @param df A provided dataframe.
#' @param metaCols Range of columns in which we have the metadata (e.g., 2:20)
#' @param colorVar A categorical variable from the dataframe to color the points.
#' (e.g., "Group"). NULL by default.
#' @importFrom dplyr select
#' @import GGally
#' @return A pairs plot.
#' @export

covPairs <- function(df, metaCols, colorVar = NULL) {
  # First, we only select the numerical variables, as the function won't
  # work with categorical variables:
  numCols <- df[, metaCols] |>
    select(where(is.numeric)) |>
    colnames()

  # We get the indexes of those columns:
  numInd <- which(colnames(df) %in% numCols)

  # We check that the colorVar variable is categorical:
  if (!is.null(colorVar) && !is.factor(df[[colorVar]])) {
    stop("colorVar must be a categorical variable.")
  }

  # And we create the plot depending on the colorVar value:
  if (is.null(colorVar)) {
    p <- ggpairs(df, columns = numInd, lower = list(continuous = "smooth"),
                 upper = list(continuous = wrap("cor", method = "spearman")))
  } else {
    p <- ggpairs(df, columns = numInd, aes(color = df[[colorVar]], alpha = 0.5),
                 lower = list(continuous = "smooth"),
                 upper = list(continuous = wrap("cor", method = "spearman")))
  }

  # Then, we print the plot:
  print(p)

  # Finally, we return the plot:
  return(p)
}
