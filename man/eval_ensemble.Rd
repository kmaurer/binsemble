% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predictFunctions.R
\name{eval_ensemble}
\alias{eval_ensemble}
\title{Function to evaluated the ensemble test accuracy against member accuracies}
\usage{
eval_ensemble(ensemble, test_data)
}
\arguments{
\item{ensemble}{An ensemble classifier built with \code{buildWeightedEnsemble}}

\item{test_data}{Data to predict with the ensemble}
}
\value{
Predictions for the provided data
}
\description{
Function to evaluated the ensemble test accuracy against member accuracies
}
