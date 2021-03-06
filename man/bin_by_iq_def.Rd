% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iqBinningFunctions.R
\name{bin_by_iq_def}
\alias{bin_by_iq_def}
\title{Iterative Quantile Binning New Data from defined bins}
\usage{
bin_by_iq_def(bin_def, new_data, output = "data", strict = FALSE)
}
\arguments{
\item{new_data}{Data frame with column names matching the binned columns from bin-training data}

\item{output}{Matches format of iterative_quant_bin and inherets properties from iqnn if applicable {"data","both"}}

\item{strict}{TRUE/FALSE: If TRUE Observations must fall within existing bins to be assigned; if FALSE the outer bins in each dimension are unbounded to allow outlying values to be assigned.}

\item{iq_def}{Iterative quantile binning definition list}
}
\value{
updated binning definition with bins extended by tolerance values
}
\description{
New observations selected from the same population as the data used to build bin definitions may fall just outside the bins. If we wish to include nearby values we can either allow outer bins to be extended (this function) or to leave the outer bins unbounded.
}
\examples{
withhold_index <- c(1,2,51,52,101,102)
iq_def <- iterative_quant_bin(data=iris[-withhold_index,], bin_cols=c("Sepal.Length","Sepal.Width","Petal.Width"),
                              nbins=c(3,2,2), output="definition")
bin_by_iq_def(bin_def=iq_def, new_data=iris[withhold_index,], output="data")
}
