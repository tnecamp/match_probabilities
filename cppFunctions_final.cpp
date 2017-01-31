#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;

/**
 * This file includes all the functions that we might use for the project. 
 * In order to source it, add the lines
 * library(Rcpp)
 * sourceCpp("name of this file")
 */

/**
 * This function is not used. 
 * Instead we use colMeans (the original one) because the cpp version did not dramatically improved timing.
 * However, we leave it here for fun.
 */

// [[Rcpp::export]]
NumericVector colMeansC(NumericMatrix x){
  NumericVector means(x.ncol());  
  int nrows = x.nrow();
  double total = 0;
  for(int i = 0; i<x.ncol(); ++i){
    
    for(int j = 0; j < nrows; ++j) {
      total += x(j,i);
    }
    means[i] = total/nrows;
    total = 0;
  }
  return means;
}

/**
 * Usage: Centers a matrix, x,  based on colMeans
 * substractColMeans(matrix to substract from, column means (via colMeans function))
 * Note that the matrix WILL be changed!
 * Example:
 * X <- matrix(runif(18, 0, 1500), ncol=3)
 * Xmean=colMeans(X) #Note that we use colMeans (the original one) because the cpp version did not dramatically improved timing.
 * substrsubstractColMeans(X, Xmean)
 */
 // [[Rcpp::export]]
void substractColMeans(NumericMatrix& x, NumericVector colMeans){
  
  int nrows = x.nrow();
  for(int i = 0; i < x.ncol(); ++i){
    
    for(int j = 0; j < nrows; ++j) {
      x(j,i) -= colMeans[i];
    }
  }
}


/**
 * Usage: Quickly computes the likelihood based on multi-variate normal
 * Inputs: likelihood_mat to store likelihood values
 * sigma_hat_inv is the precision matrix
 * test_means is the means for the multi-variate normal distirbution 
 * (test_means is based on OLS estimation and covariates)
 * Y is our data for which to compute the likelihood
 * n_test is the number of individuals in survey data
 * n_train is the number of individuals big_data
 * VB_dim is the number of covariates in big_data
 */
 // [[Rcpp::depends(RcppEigen)]]
 using Eigen::MatrixXd;
 using Eigen::VectorXd;
 // [[Rcpp::export]]
NumericMatrix likelihoodComp(NumericMatrix& likelihood_mat, 
  MatrixXd& sigma_hat_inv, MatrixXd& test_means, NumericMatrix& centered_v_b,
  int n_test, int n_train, int VB_dim){
  MatrixXd cur_c(1,VB_dim);
  MatrixXd cur_b(1, VB_dim);
  MatrixXd diff(1, VB_dim);
  double val;
  for(int i = 0; i < n_test; ++i){
    for(int k = 0; k < VB_dim; ++k){
      cur_c(0,k) = test_means(i,k);
    }
    for(int j = 0; j < n_train; ++j){
        for(int b = 0; b < VB_dim; ++b){
        cur_b(0,b) = centered_v_b(j,b);
      }
      diff = cur_b - cur_c;
      val = (-.5*diff*sigma_hat_inv*diff.transpose())(0,0);  //change to give you float
      likelihood_mat(i,j) = val;
    }
  }
  return likelihood_mat;
}


/**
 * Usage: Quickly computes eigen vector and eigenvalues of a matrix
 */
#include <Eigen/Eigenvalues>
 using Eigen::MatrixXcd;
  // [[Rcpp::depends(RcppEigen)]]
 // [[Rcpp::export]]
MatrixXcd eigenSolving(MatrixXd Q){
  Eigen::EigenSolver<MatrixXd> es(Q);
  MatrixXcd A = es.eigenvectors();
  MatrixXcd B = es.eigenvalues();
  MatrixXcd C(A.rows()+B.rows(), A.cols());
  C << A,
       B;
  return C;
}


