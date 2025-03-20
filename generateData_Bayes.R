################################################################################
############################## generateData_Bayes.R ############################
################################################################################

# Inputs: As inputs it takes the output values of Model.R
#         mu1, mu2, sigma1, sigma2, Omega1, Omega2, logdet and
#         dimension (p), train/test sample sizes of population 1: n1/m1, and
#         population 2: n2/m2

#         It then generates the data from multivariate gaussian distribution 
#         using the specific mu1, mu2, sigma1 and sigma2 and calculates the  
#         class assignments using the bayes rule.
# 
# Outputs: Data: X1 ((n1 + m1) x p dim matrix), 
#                X2 ((n2 + m2) x p dim matrix)
#          res: Estimated sequence of 1s and 2s using bayes rule. 


################################################################################
################################################################################

generateData_Bayes <- function(mu1, mu2, sigma1, sigma2, Omega1, Omega2, logdet, n1 = 100, n2 = 100, m1 = 100, m2 = 100) {
  
  library(MASS)
  library(Matrix)
  library(pracma)
  
  p <- length(mu1)
    
    X1 <- as.matrix(t(mu1 + crossprod((chol(sigma1)), matrix(rnorm((n1 + m1)*p),nrow = p,ncol = (n1 + m1)))))
    X2 <- as.matrix(t(mu2 + crossprod((chol(sigma2)), matrix(rnorm((n2 + m2)*p),nrow = p,ncol = (n2 + m2)))))
    
    term1 <- c(0)
    term2 <- c(0)
    X_test <- rbind(X1[(n1 + 1):(n1 + m1),], X2[(n2 + 1):(n2 + m2),])
    for(i in 1:(m1 + m2))
    {
      term1[i] <- tcrossprod((tcrossprod((X_test[i,]-t(mu1)),Omega1)),((X_test[i,]-t(mu1))))
      term2[i] <- tcrossprod((tcrossprod((X_test[i,]-t(mu2)),Omega2)),((X_test[i,]-t(mu2))))
    }
    dis <- term2 - term1 + logdet + 2*log(m1/m2)
    res <- ifelse(dis>0, 1, 2)
    
    return(list(X1 = X1, X2 = X2, res = res ))
}
