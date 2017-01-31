\name{matchgraph}
\alias{matchgraph}
\title{
Generate histogram of the output from match_prabability function
}
\description{
This function returns the histogram of matching probabilities of a certain individual.
}
\usage{
matchgraph(probability_matrix, survey_ind, min_prob = NULL, 
upper_quantile = NULL)
}
\arguments{
  \item{probability_matrix}{
  Matching probability table, usually from the output of match_probability function.
}
  \item{survey_ind}{
Index of row number of the desired individual.
}
  \item{min_prob}{
  All probabilities that are larger than it will be counted. If min_prob=NULL, all probabilities will be counted.
}
  \item{upper_quantile}{
  Value must be between 0 and 1. if assigned, the probabilities above the upper quantile will be counted.
}
}
\details{
Note that min_prob and upper_quantile cannot be assigned values at the same time.
}
\value{
This function will return a histogram of the matching probabilities of the assigned individual with conditions given by min_prob or upper_quantile.
}
\author{
Timothy NeCamp, Dana Turjeman, Xiaomeng Du
}
\examples{
##create example data
a<-matrix(seq(0.001,0.1,by=0.001),nrow=4)

matchgraph(a,1,min_prob=0.01)
matchgraph(a,1,upper_quantile=0.5)

}