\name{auto_poly_koyck}
\alias{auto_poly_koyck}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
auto_poly_koyck: Automatic Dynamically Complete Koyck Lag Model
}
\description{
Takes a dependent variable and potentially multiple independent variables to construct a dynamically complete koyck distributed lag model. Also tests for the stability of the model's lag function.
}
\usage{
auto_poly_koyck(data, dependent_variable, independent_variables, max_lags = 10, bg_test_level = 0.10, L = "L", print_output = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{---Data frame that variable inputs reside within.}
  \item{dependent_variable}{---The numerical variable which will be used as the dependent variable and of which lags will be created and added into the linear model.}
  \item{independent_variables}{---A numerical vector composed of all of the independent variables. See "examples" for how they are to be inputed.}
  \item{max_lags}{---Numerical argument outlining the maximum number of lagged dependent variable that can be included in the linear model. Default is 10.}
  \item{L}{---Parameter used for polynomial function creation.}
  \item{bg_test_level}{---Numerical argument outlining the significance level that the Breusch-Godfrey test must exclipse before termination of iterative lag addition. Default is 0.10.}
  \item{print_output}{---Logical argument determining if printing of polynomial root testing is printed.}

}
\details{
Model selection is done via iterative addition of lagged dependent variables until the Breuch-Godfrey test of order 1 yields a p-value greater than 0.10: removes residual serial correlation.
The resulting parameter coefficients of the lagged dependent variables are then used to test if the resulting polynomial lag function has roots inside the unit circle, which would result in an "explosive" lag function.
Model uses lm() and can handle missing values automatically through listwise deletion.
}

\value{
  \item{real}{---list comprising of the real components of the derived polynomial roots: numeric.}
  \item{complex}{---list comprising of the imaginary components of the derived polynomial roots: numeric.}
  \item{outside}{---list denoting the unit roots status of each polynomial root: logical.}

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Solomon C. Salama - University of London
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
'bgtest()'
}
\examples{


{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
