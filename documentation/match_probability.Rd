\name{match_probability}
\alias{match_probability}
\title{
Matching probability
}
\description{
This main function calculates the matching probability table of the two datasets.
}
\usage{
match_probability(big_data, survey_data, common_var_index)
}
\arguments{
  \item{big_data}{
A big dataset that contains two kinds of variables in it: those are also in survey_data, and those that are unique. Must can be trasformed into Matrix.
}
  \item{survey_data}{
A smaller dataset that contains common variables that are also in big_data. Could contain unique variables that won't be used. Must can be trasformed into Matrix.
}
  \item{common_var_index}{
An n by 2 matrix that contains the matching index of common variables of the two datasets. The first column stands for the index of common variables in big_data, while the second column in survey_data. The number of columns will always be two, while the number of rows equals to the total number of common variables.
}
}
\value{
The output of the main function is a matching probability table with the number of columns equals to the number of rows in big_data, and number of rows equals to the number of rows in survey_data. Each row of this table contains the probability of the individual in survey_data being matched with the individuals in big_data.
}
\author{
Timothy NeCamp, Dana Turjeman, Xiaomeng Du
}
\examples{
##creating example inputs
n_train <- 100
n_test <- 5
VC_dim <- 5
VB_dim <- 2
sigma_VC <- 1
sigma_noise <- .5
VB_corr <- .2

VC<-matrix(rnorm((n_train+n_test)*VC_dim,0,sigma_VC),ncol=VC_dim)
beta <- matrix(rnorm(VC_dim*VB_dim, sd = 10), nrow = VC_dim, ncol = VB_dim)
##beta<-matrix(rep(1,VC_dim*VB_dim),nrow=VC_dim, ncol = VB_dim)
VC_test <- VC[1:n_test,]
VC <- VC[(n_test+1):(dim(VC)[1]),]
dim(VC_test) <- c(n_test, VC_dim)
dim(VC) <- c(n_train, VC_dim)
VC_test[1:2,] <- VC[1:2,]


sigma<-matrix(VB_corr, ncol=VB_dim, nrow=VB_dim)
diag(sigma) <- 1
sigma <- sigma*sigma_noise
VB<-matrix(nrow=n_train,ncol=VB_dim)


temp <- VC\%*\%beta
for (i in 1:dim(VC)[1]){
VB[i,]<-mvrnorm(mu=(temp)[i,],Sigma=sigma)}
VB <- as.matrix(VB)

big_data <- cbind(VB,VC)
survey_data <- cbind(VC_test, matrix(rep(1), nrow = n_test, ncol=3))
common_var_index <- cbind((VB_dim+1):(VC_dim+VB_dim),1:VC_dim)

##using match_probability to calculate matching table
match_prob <- match_probability(big_data, survey_data, common_var_index)
}
