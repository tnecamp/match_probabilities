##This main function calculates the match probability table of the two datasets ##

match_probability = function(big_data, survey_data, common_var_index){
	#big_data -A big dataset that contains two kinds of variables in it: those are also in survey_data, and those that are unique
    #survey_data- A smaller dataset that contains common variables that are also in big_data. Could contain unique variables that won't be used. Must can be trasformed into Matrix
	#common_var_index -An n by 2 matrix that contains the matching index of common variables of the two datasets. The first column stands for the index of common variables in big_data, while the second column in survey_data. The number of columns will always be two, while the number of rows equals to the total number of common variables.	
	
	n_test = dim(survey_data)[1]	# length of survey data
	n_train = dim(big_data)[1]	# length of big_data
	big_data_common = common_var_index[,1]	# index of common variables in big data
	survey_common = common_var_index[,2]	# index of common variables in survey data
	VC = big_data[,big_data_common]		# matrix of big data common variables
	VC_test = survey_data[,survey_common]	 # matrix of survey data common variable
	VB = big_data[,-big_data_common]	# matrix of variables in bgi data not in common variables
	
	# turn everything into matrices
	X=as.matrix(VC)
	Y=as.matrix(VB)
	VC_test = as.matrix(VC_test)
	
	VB_dim = dim(Y)[2]
	VC_dim = dim(X)[2]
	
	# center our data
	Xmean=colMeans(X)
	substractColMeans(X, Xmean)
	Ymean=colMeans(Y)
	substractColMeans(Y, Ymean)
	
	# get matrix which shares information across OLS regression on each dimension based on 
	# eigen values, see write-up for further details
	if(VB_dim > 1){
	M=t(Y)%*%X
	Q=solve(t(Y)%*%Y)%*%M%*%solve(t(X)%*%X)%*%t(M)
	eigen_mat = eigenSolving(Q)
	eigen_mat = Re(eigen_mat)
	C2 = eigen_mat[(VB_dim+1):(2*VB_dim) , 1]
	T_mat = eigen_mat[1:VB_dim, 1:VB_dim]
	T_mat = t(T_mat)
	
	r=dim(X)[2]/dim(X)[1]
	D=(1-r)*(C2-r)/((1-r)^2*C2+r^2*(1-C2))
	D=diag(replace(D,D<0,0))
	B=solve(T_mat)%*%D%*%T_mat
	} else {
		B = 1
		}
	
	# get OLS estimates for each column of big_data
	theta1=matrix(NA,nrow=dim(X)[2],ncol=dim(Y)[2])
	for (j in 1:(dim(Y)[2])){
		theta1[,j]=fastLm(Y[,j]~0+X)$coefficients  #obtain OLS estimates using faster algorithm
	}

	
	centered_v_b = Y
	centered_v_c = X
	theta_hat = theta1
	b_star = B
	mean_v_c = Xmean
	v_c_test = VC_test
	
	# compute estimate of covariance matrix for likelihood
	compute_sigma_hat = function(centered_v_b, centered_v_c, theta_hat, b_star){
		resid = centered_v_b - centered_v_c%*%theta_hat%*%b_star
		return(cov(resid))
		}
	
	sigma_hat = compute_sigma_hat(Y, X, theta1, B)
	
	# compute the likelhiood based on likelihood of interest. This is based on multi-variate
	# normal likelihood
	substractColMeans(VC_test, Xmean)
	test_means = VC_test%*%theta1%*%B
	sigma_hat_inv = solve(sigma_hat)
	likelihood_mat_c = matrix(nrow=n_test, ncol = n_train)

	likelihood_mat = likelihoodComp(likelihood_mat_c, sigma_hat_inv, test_means, Y, n_test, n_train, VB_dim)		# function in C which quickly computes multi-variate normal log likelihood
	likelihood_mat = exp(likelihood_mat)
	
	# normalize to get probabilities
	denom = rowSums(likelihood_mat)
	denom[denom==0] = 10
	prob_mat = likelihood_mat/denom  	
	
	# return the matrix of probabilities of matching, where entry i,j is the probability of 
	# matching individual i in the survey_data to individual j in the big_data
	return(prob_mat) 
	}
	

## topmatch is a function to get the top n matches of each individual  ##
## in the survey data from a given probability matrix ##

topmatch=function(probability_matrix, survey_ind, n=5){
	# probability matrix is the outpout from match_probability
	# survey_ind is the individual for whom you want the histrogram
	# n specifies the number of top matches you wish to include
	
	x = survey_ind
	y=n
	probability_matrix[x,order(probability_matrix[x,],decreasing=T)[1:y]]
	}


## Generates histogram of the output from match_prabability for a certain survey individual ##

matchgraph=function(probability_matrix, survey_ind, min_prob=NULL, upper_quantile=NULL){
	# probability matrix is the outpout from match_probability
	# survey_ind is the individual for whom you want the histrogram
	# min_prob is the minimal probability you want included in your histogram
	# upper_quantile forces the histogram to only include probabilities above the upper quantile will be counted

	x = survey_ind
	p = min_prob
	alpha = upper_quantile
	if(is.null(p)&&is.null(alpha)){
		hist(probability_matrix[x,],main=paste("Matching Probability of Individual ",x),xlab="matching probability")
		}else if(!is.null(p)&&is.null(alpha)){
			hist(probability_matrix[x,probability_matrix[x,]>p],main=paste("Matching Probability of Individual ",x," Greater Than ",p),xlab="matching probability")
			}else if(is.null(p)&&!is.null(alpha)){
				hist(probability_matrix[x,probability_matrix[x,]>quantile(probability_matrix[x,],(1-alpha))],main=paste("Matching Probability of Individual ",x," Top ",alpha*100,"%"),xlab="matching probability")
				}else print("Error: p and alpha cannot be NULL at the same time")
	}


## This function shows the overall matching probabilities (aggregate over all users)  ##
## to be used  to make decision about whether it is a good chance of match. ##
overallmatch=function(probability_matrix, min_prob = NULL,n = 10){
	# probability matrix is the outpout from match_probability
	# min_prob is the minimal probability you want included in your histogram
	# n is the number of top probabilities
	
	p = min_prob
	if(is.null(p)){
	return(probability_matrix[order(probability_matrix,decreasing=T)[1:n]])
		}else{
		return(probability_matrix[probability_matrix[]>p])
		}
	}