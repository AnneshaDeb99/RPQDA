Bayes <- function(X_test,mu1,mu2,sigma1,sigma2,n1,n2){

  term3 <- 2*sum(log(diag(chol(sigma2))))- 2*sum(log(diag(chol(sigma1))))
  m <- nrow(X_test)
  term1 <- c(0)
  term2 <- c(0)
  for(i in 1:m)
  {
    term1[i] <- (X_test[i,]-t(mu1))%*%chol2inv(chol(sigma1))%*%(t(X_test[i,]-t(mu1)))
    term2[i] <- (X_test[i,]-t(mu2))%*%chol2inv(chol(sigma2))%*%(t(X_test[i,]-t(mu2)))
  }
  dis <- term2-term1+term3 + log(n1/n2)
  res <- ifelse(dis>0,0,1)
  return(res)
}

# logdet <- function(A) {
#   ans <- 2*sum(log(diag(chol(A))))
#   return(ans)
# }
