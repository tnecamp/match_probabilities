\name{topmatch}
\alias{topmatch}
\title{
Generate top probability of one individual
}
\description{
This function returns the top n matching probability of the survey_ind individual in Survey_Data.
}
\usage{
topmatch(probability_matrix, survey_ind, n = 5)
}
\arguments{
  \item{probability_matrix}{
Matching probability table.
}
  \item{survey_ind}{
Index of row number of the desired individual.
}
  \item{n}{
The number of top probabilities, by default equals to 5.
}
}
\value{
The output is a vector of the top n matching probabilities for survey respondent that is assigned.
}
\author{
Xiaomeng Du
}
\examples{

a<-matrix(seq(0.001,0.1,by=0.001),nrow=4)
topmatch(a,1)

}