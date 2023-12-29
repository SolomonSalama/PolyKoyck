\name{adf_kpss_test}
\alias{adf_kpss_test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
adf_kpss_test: Joint Testing for Stationarity and Classification
}
\description{
Uses Augmented Dickey-Fuller (tseries) and KPSS test (tseries) to determine the stationarity of a time series variable, and plots PACF of variable.
Lag length of ADF test is determined by the number of significant PACF lags before an insignificant one is detected.
}
\usage{
adf_kpss_test(time_series_variable)
%- maybe also 'usage' for other objects documented here.
}
\arguments{
  \item{time_series_variable}{---Numeric variable with no missing values.}
%%     ~~Describe \code{x} here~~
}

\details{
This function sequentially goes through the lags of the pacf plot from the pacf() function. Beginning with lag one, it adds up the lags which are statistically significant (95 percent CI) until it finds an insignificant one. If no lags beyond lag one are significant, the adf.test() takes k = 1. The KPSS test set the null to "Trend" and lshort = TRUE. There are four potential results: stationary, trend-stationary, difference-stationary, and non-stationary. Does not take missing values.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{

Kwiatkowski, D.; Phillips, P. C. B.; Schmidt, P.; Shin, Y. (1992). "Testing the null hypothesis of stationarity against the alternative of a unit root". Journal of Econometrics. 54 (1–3): 159–178.

S. E. Said and D. A. Dickey (1984). "Testing for Unit Roots in Autoregressive-Moving Average Models of Unknown Order". Biometrika. 71: 599-607.

https://www.statsmodels.org/dev/examples/notebooks/generated/stationarity_detrending_adf_kpss.html

https://www.rdocumentation.org/packages/tseries/versions/0.10-55/topics/adf.test

}
\author{
Solomon C. Salama - University of London
}
\note{
Example code pulled from R documentation page of adf.test() from the 'tseries' library
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
'pacf()'
'adf.test()'
'kpss.test()'
}
\examples{
x <- rnorm(100)  #stationarity process
adf_kpss_test(x)

y <- cumsum(x)  #non-stationary process
adf_kpss_test(y)

x <- 0.8*(1:100)+rnorm(100)  #stationarity process
adf_kpss_test(x)
}


