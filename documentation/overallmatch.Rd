\name{overallmatch}
\alias{overallmatch}
\title{
Overall matching quality
}
\description{
This function shows the overall matching probabilities and users can use it to make decision about whether it is a good chance of match.
}
\usage{
overallmatch(probability_matrix, min_prob = NULL, n = 10)
}
\arguments{
  \item{probability_matrix}{
  Matching probability table.
}
  \item{min_prob}{
  All probabilities that are larger than it will be counted. If min_prob=NULL, all probabilities will be counted.
}
  \item{n}{
  The number of top probabilities, by default equals to 10.
}
}
\value{
This output shows the top n matching probabilities regardless of which respondent we calculate. It can tell us the best matches of the two datasets. 
}
\author{
Timothy NeCamp, Dana Turjeman, Xiaomeng Du
}
\examples{
a<-matrix(seq(0.001,0.1,by=0.001),nrow=4)
overallmatch(a)

}