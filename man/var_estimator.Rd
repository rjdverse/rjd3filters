% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/crossvalidation.R
\name{var_estimator}
\alias{var_estimator}
\title{Variance Estimator}
\usage{
var_estimator(x, coef, ...)
}
\arguments{
\item{x}{input time series.}

\item{coef}{vector of coefficients or a moving-average (\code{\link[=moving_average]{moving_average()}}).}

\item{...}{other arguments passed to the function \code{\link[=moving_average]{moving_average()}} to convert \code{coef} to a \code{"moving_average"} object.}
}
\description{
Variance Estimator
}
\details{
Let \eqn{(\theta_i)_{-p\leq i \leq q}} be a moving average of length \eqn{p+q+1} used
to filter a time series \eqn{(y_i)_{1\leq i \leq n}}.
It is equivalent to a local regression and the associated error variance \eqn{\sigma^2}
can be estimated using the normalized residual sum of squares, which can be simplified as:
\deqn{
\hat\sigma^2=\frac{1}{n-(p+q)}\sum_{t=p+1}^{n-q}
\frac{(y_t-\hat \mu_t)^2}{1-2w_0^2+\sum_{i=-p}^q w_i^2}
}
}
\references{
Loader, Clive. 1999.
Local regression and likelihood.
New York: Springer-Verlag.
}
\seealso{
\code{\link[=df_var]{df_var()}}.
}
