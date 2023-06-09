% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/covAnalysis.R
\name{covAnalysis}
\alias{covAnalysis}
\title{covAnalysis}
\usage{
covAnalysis(
  df,
  metaCols,
  patientCol,
  bioCols,
  ind_cov,
  REML,
  printFormula = TRUE,
  resultsFolder = "results"
)
}
\arguments{
\item{df}{A given dataframe.}

\item{metaCols}{Range of columns with metadata information (e.g., 2:15)}

\item{patientCol}{Name of the ID or patient column in the top-table.}

\item{bioCols}{Range of columns with biomarker information (e.g., 16:300)}

\item{ind_cov}{Individual covariates for which we want to test the interaction
with the biomarkers. They will be included as fixed effects in all models.}

\item{REML}{Algorithm used for estimating the model parameters. It is recommended
to set it as FALSE (it will use ML) if we are comparing models.}

\item{printFormula}{Print the formula of the different models (TRUE by default).}

\item{resultsFolder}{Name of the directory in which the results will be stored
("results" by default).}
}
\value{
A dataframe with the results of the mediation analysis
}
\description{
This function analyzes the influence of the covariates of a given dataframe
in the levels of the selected biomarkers using linear mixed effects models.
The main purpose of this function is to detect the existence of potential
confounding covariates. It takes several arguments as input: a given
dataframe ("df"), the range of columns with the covariates (metaCols),
the name of the patient or ID column (patientCol), the range of columns
that contain the biomarker data (bioCols), the individual covariates which
will be included as an interaction the model (ind_cov), the method used to
evaluate the models (REML) and an optional argument (TRUE by default) that
will print the formula of each model.
}
\details{
The function iterates through each metadata variable (potential confounding
variables), generates the formulas for the linear mixed-effects models, and
calls the lmeBioCov function that will generate a model for all biomarkers.
The results from each iteration are added to the results dataframe. Finally,
the results are exported to an Excel file and the function returns the
a dataframe with all the models generated.
}
