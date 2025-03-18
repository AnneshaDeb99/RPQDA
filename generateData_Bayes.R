################################################################################
############################## generateData_Bayes.R ############################
################################################################################

# Inputs: As inputs it takes the output values of Model.R
#         mu0, mu1, sigma0, sigma1, Omega0, Omega1, logdet and
#         dimension (p), train/test sample sizes of population 1: n0/m0, and
#         population 2: n1/m1
#         It then generates the data from multivariate gaussian distribution 
#         using the specific mu0, mu1, sigma0 and sigma1 and calculates the  
#         class assignments using the bayes rule.
# 
# Outputs: Data: X0 ((n0+m0)x p dim matrix), 
#                X1 ((n1+m1)x p dim matrix)
#          res: Estimated sequence of 0s and 1s using bayes rule. 


################################################################################
################################################################################

generateData_Bayes <- function(mu0, mu1, sigma0, sigma1, Omega0, Omega1, logdet , p, n0 = 100, n1 = 100,m0 = 100, m1 = 100) {
  
  #Working
  
  library(MASS)
  library(Matrix)
  library(pracma)
    
    X0 <- as.matrix(t(mu0 + crossprod((chol(sigma0)), matrix(rnorm((n0+m0)*p),nrow=p,ncol=(n0+m0)))))
    X1 <- as.matrix(t(mu1 + crossprod((chol(sigma1)), matrix(rnorm((n1+m1)*p),nrow=p,ncol=(n1+m1)))))
    
    term0 <- c(0)
    term1 <- c(0)
    X_test <- rbind(X0[(n0+1):(n0+m0),],X1[(n1+1):(n1+m1),])
    for(i in 1:(m0+m1))
    {
      term0[i] <- tcrossprod((tcrossprod((X_test[i,]-t(mu0)),Omega0)),((X_test[i,]-t(mu0))))
      term1[i] <- tcrossprod((tcrossprod((X_test[i,]-t(mu1)),Omega1)),((X_test[i,]-t(mu1))))
    }
    dis <- term1-term0+ logdet + 2*log(m0/m1)
    res <- ifelse(dis>0,0,1)
    
    return(list(X0 = X0, X1 = X1,res = res ))
}
