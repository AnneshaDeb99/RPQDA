################################################################################
################################# RPE.R ########################################
################################################################################

# Inputs: type: distribution of random matrix to be used (e.g., 'StdNormal', '3_point' etc)
#         X_train: n x p dimensional data matrix (covariates of training observations)
#         y_train: classes corresponding to the training data
#         X_test: m x p dimensional data matrix (covariates of test observations)
#         d: reduced dimension
#         B: number of random matrices to be used
#
# Outputs: class.RPE: estimated classes 

################################################################################
################################################################################
################################################################################

library(abind)

RPE <- function(type, X_train, y_train, X_test, d = 10, B = 500){
  
  if(is.vector(X_test) == TRUE){
    X_test = t(as.matrix(X_test))
  }
  k <- length(unique(y_train)) #number of classes 
  p <- ncol(X_test)
  n.test <- nrow(X_test)
  
  if(type == 'StdNormal'){
    q_RP_hat_ens <- foreach(b = 1:B, .combine = function(...) abind(..., along = 3), .multicombine = TRUE, .init = NULL, .packages = c('doParallel','abind','Matrix')) %dopar% {
      
      set.seed(2024 + b)
      R <- matrix(rnorm(d*p,0,1), nrow = p, ncol = d)
      XRP_train <- crossprod(t(X_train), R)
      XRP_test <- crossprod(t(X_test), R)
      n.train <- nrow(XRP_train)
      n.test <- nrow(XRP_test)
      q.RP.hat = matrix(0, ncol = k, nrow = n.test)
      
      for(i in 1:n.test){
        
        xrp.te = XRP_test[i, ]
        q.RP.hat[i, ] = sapply(1:k, function(j){
          n = sum(y_train == j);
          class = which(y_train == j);
          xrp.train = XRP_train[class, ];
          mu = colMeans(xrp.train);
          sig = crossprod(xrp.train - mu, xrp.train - mu)/(n-1);
          dis = - (0.5*log(abs(det(sig)))) - (0.5*tcrossprod((tcrossprod(xrp.te - mu, solve(sig))), xrp.te - mu)) + (log(n/n.train));
          return(dis)})
      }
      return(q.RP.hat)
    }
  } 
  else if(type == 'Normal'){
    q_RP_hat_ens <- foreach(b = 1:B, .combine = function(...) abind(..., along = 3), .multicombine = TRUE, .init = NULL, .packages = c('doParallel','abind')) %dopar% {
      
      set.seed(2024 + b)
      R <- matrix(rnorm(d*p,0,(1/d)),nrow=p,ncol=d)
      XRP_train <- crossprod(t(X_train), R)
      XRP_test <- crossprod(t(X_test), R)
      n.train <- nrow(XRP_train)
      n.test <- nrow(XRP_test)
      q.RP.hat = matrix(0, ncol = k, nrow = n.test)
      
      for(i in 1:n.test){
        
        xrp.te = XRP_test[i, ]
        q.RP.hat[i, ] = sapply(1:k, function(j){
          n = sum(y_train == j);
          class = which(y_train == j);
          xrp.train = XRP_train[class, ];
          mu = colMeans(xrp.train);
          sig = crossprod(xrp.train - mu, xrp.train - mu)/(n-1);
          dis = - (0.5*log(det(sig))) - (0.5*tcrossprod((tcrossprod(xrp.te - mu, solve(sig))), xrp.te - mu)) + (log(n/n.train));
          return(dis)})
      }
      return(q.RP.hat)
    }
  }else if(type == 'Std+1/-1'){
    q_RP_hat_ens <- foreach(b = 1:B, .combine = function(...) abind(..., along = 3), .multicombine = TRUE, .init = NULL, .packages = c('doParallel','abind')) %dopar% {
      
      set.seed(2024 + b)
      r <- rbinom(d*p, 1, 0.5)
      r[r == 0] <- -1
      R <- matrix(r/d , nrow = p, ncol = d)
      
      XRP_train <- crossprod(t(X_train), R)
      XRP_test <- crossprod(t(X_test), R)
      n.train <- nrow(XRP_train)
      n.test <- nrow(XRP_test)
      q.RP.hat = matrix(0, ncol = k, nrow = n.test)
      
      for(i in 1:n.test){
        
        xrp.te = XRP_test[i, ]
        q.RP.hat[i, ] = sapply(1:k, function(j){
          n = sum(y_train == j);
          class = which(y_train == j);
          xrp.train = XRP_train[class, ];
          mu = colMeans(xrp.train);
          sig = crossprod(xrp.train - mu, xrp.train - mu)/(n-1);
          dis = - (0.5*log(det(sig))) - (0.5*tcrossprod((tcrossprod(xrp.te - mu, solve(sig))), xrp.te - mu)) + (log(n/n.train));
          return(dis)})
      }
      return(q.RP.hat)
    }
  }
  else if(type == '+1/-1'){
    q_RP_hat_ens <- foreach(b = 1:B, .combine = function(...) abind(..., along = 3), .multicombine = TRUE, .init = NULL, .packages = c('doParallel','abind')) %dopar% {
      
      set.seed(2024 + b)
      r <- rbinom(d*p, 1, 0.5)
      r[r == 0] <- -1
      R <- matrix(r, nrow = p, ncol = d)
      
      XRP_train <- crossprod(t(X_train), R)
      XRP_test <- crossprod(t(X_test), R)
      n.train <- nrow(XRP_train)
      n.test <- nrow(XRP_test)
      q.RP.hat = matrix(0, ncol = k, nrow = n.test)
      
      for(i in 1:n.test){
        
        xrp.te = XRP_test[i, ]
        q.RP.hat[i, ] = sapply(1:k, function(j){
          n = sum(y_train == j);
          class = which(y_train == j);
          xrp.train = XRP_train[class, ];
          mu = colMeans(xrp.train);
          sig = crossprod(xrp.train - mu, xrp.train - mu)/(n-1);
          dis = - (0.5*log(det(sig))) - (0.5*tcrossprod((tcrossprod(xrp.te - mu, solve(sig))), xrp.te - mu)) + (log(n/n.train));
          return(dis)})
      }
      return(q.RP.hat)
    }
  }  
  else if(type == '3_point'){
    q_RP_hat_ens <- foreach(b = 1:B, .combine = function(...) abind(..., along = 3), .multicombine = TRUE, .init = NULL, .packages = c('doParallel','abind')) %dopar% {
      
      set.seed(2024 + b)
      R <- matrix(sample(c(1,0,-1),(d*p),replace = TRUE, prob = c(1/(2*sqrt(p)), 1 - 1/sqrt(p), 1/(2*sqrt(p)))), nrow = p, ncol = d)
      
      XRP_train <- as.matrix(crossprod(t(X_train), R))
      XRP_test <- as.matrix(crossprod(t(X_test), R))
      n.train <- nrow(XRP_train)
      n.test <- nrow(XRP_test)
      q.RP.hat = matrix(0, ncol = k, nrow = n.test)
      
      for(i in 1:n.test){
        
        xrp.te = XRP_test[i, ]
        q.RP.hat[i, ] = sapply(1:k, function(j){
          n = sum(y_train == j);
          class = which(y_train == j);
          xrp.train = XRP_train[class, ];
          mu = colMeans(xrp.train);
          sig = crossprod(xrp.train - mu, xrp.train - mu)/(n-1);
          dis = - (0.5*log(abs(det(sig)))) - (0.5*tcrossprod((tcrossprod(xrp.te - mu, solve(sig))), xrp.te - mu)) + (log(n/n.train));
          return(dis)})
      }
      return(q.RP.hat)
    }
  }
  class.RPE = sapply(1:n.test, function(n){
    discriminant = rowMeans(q_RP_hat_ens[n,,]);
    class = which(discriminant == max(discriminant))
    return(class)
  })
  return(class.RPE)
}
