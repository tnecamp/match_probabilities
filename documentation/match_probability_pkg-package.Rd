\name{matchProbabilityPkg-package}
\alias{matchProbabilityPkg-package}
\alias{match_probability_pkg}
\docType{package}
\title{
  Calculating matching probabilities of two datasets
}
\description{
  matchProbabilityPkg gives a clear measure of what are the probabilities that any one individual in one database can be identified as being any specific individual from another database. 
}
\details{
  Our package contains four functions. The main function calculates the matching probability table, and the rest of the functions use this table to generate summaries of the table.
  Main function: match_probability()
  Other functions:  topmatch(), matchgraph(),overallmatch().
  For more details, please check the help pages of each function.
}
\author{
Timothy NeCamp, Dana Turjeman, Xiaomeng Du

Maintainer: Timothy NeCamp, Dana Turjeman, Xiaomeng Du
}
\references{
  Breiman, Leo, and Jerome H. Friedman. "Predicting multivariate responses in multiple linear regression." Journal of the Royal Statistical Society: Series B (Statistical Methodology) 59.1 (1997): 3-54.
  
  Gilula, Z., McCulloch, R. E., & Rossi, P. E. (2006). A direct approach to data fusion. Journal of Marketing Research, 43(1), 73-83. 
  
  Ong, A. D., & Weiss, D. J. (2000). The impact of anonymity on responses to sensitive Questions1. Journal of Applied Social Psychology, 30(8), 1691-1708. 
  
Qian, Y., & Xie, H. (2013). Which brand purchasers are lost to counterfeiters? an application of new data fusion approaches. Marketing Science, 33(3), 437-448. 
}