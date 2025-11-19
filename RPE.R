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
#Note: RPE-SN and RPE-TP are RPE with type: 'StdNormal' and '3_point' respectively.
################################################################################
################################################################################
################################################################################

logdet = function(A) {
  L = eigen(A)$values
  a = sum(log(L))
  return(a)
}

discriminant = function(X_test, mu1_hat, mu2_hat, sig1_hat, sig2_hat, n1, n2){
  sig1_inv = solve(sig1_hat)
  sig2_inv = solve(sig2_hat)
  term3 = logdet(sig2_hat)-logdet(sig1_hat)
  m = nrow(X_test)
  term1 = term2 = c(0)
  for(i in 1:m)
  {
    term1[i] = tcrossprod(tcrossprod((X_test[i,]-t(mu1_hat)),sig1_inv),(X_test[i,]-t(mu1_hat)))
    term2[i] = tcrossprod(tcrossprod((X_test[i,]-t(mu2_hat)),sig2_inv),(X_test[i,]-t(mu2_hat)))
  }
  dis = term2 - term1 + term3 + log(n1/n2)
  return(dis)
}



RPE <- function(type, X_train, y_train, X_test, d = 10, B = 500){
  ##Data
  if(is.matrix(X_test) == FALSE){
    X_test = (as.matrix(X_test))
  }
  train = cbind(y_train,X_train)
  train1 = train[train[,1]==1,]
  train2 = train[train[,1]==2,]
  n1 = length(train1[,1])   #train obs for class 1
  n2 = length(train2[,1])            #train obs for class 2
  m = nrow(as.matrix((X_test)))   #test observations
  p = ncol(as.matrix((X_train)))  #original dimension

  if(type == 'StdNormal'){
    q_RP_hat1 = foreach(b = 1:B, .combine = 'cbind', .export = c('discriminant','logdet'), .packages = 'doParallel') %dopar% {
      R = matrix(rnorm(d*p,0,1),nrow=p,ncol=d)
      XRP_train = crossprod(t(X_train),R)
      train = cbind(y_train,XRP_train)
      train1 = train[train[,1]==1,]
      train2 = train[train[,1]==2,]
      mu1 = colMeans(train1[,-1])
      mu2 = colMeans(train2[,-1])
      sig1 = crossprod((train1[,-1]-mu1),(train1[,-1]-mu1))/(n1-1)
      sig2 = crossprod((train2[,-1]-mu2),(train2[,-1]-mu2))/(n2-1)
      q_RP_hat = discriminant(X_test%*%R,mu1,mu2,sig1,sig2,n1,n2)
      return(q_RP_hat)
    }
  }
  else if(type == 'Normal'){
    q_RP_hat1 = foreach(b = 1:B, .combine = 'cbind', .export = c('discriminant','logdet')) %dopar% {
      R = matrix(rnorm(d*p,0,(1/d)),nrow=p,ncol=d)
      XRP_train = crossprod(t(X_train),R)
      train = cbind(y_train,XRP_train)
      train1 = train[train[,1]==1,]
      train2 = train[train[,1]==2,]
      mu1 = colMeans(train1[,-1])
      mu2 = colMeans(train2[,-1])
      sig1 = crossprod((train1[,-1]-mu1),(train1[,-1]-mu1))/(n1-1)
      sig2 = crossprod((train2[,-1]-mu2),(train2[,-1]-mu2))/(n2-1)
      q_RP_hat = discriminant(X_test%*%R,mu1,mu2,sig1,sig2,n1,n2)
      return(q_RP_hat)
    }
  }
  else if(type == 'Std+1/-1'){
    q_RP_hat1 = foreach(b = 1:B, .combine = 'cbind', .export = c('discriminant','logdet')) %dopar% {
      r = rbinom(d*p,1,0.5)
      r[r==0] = -1
      R = matrix(r,nrow=p,ncol=d)
      XRP_train = crossprod(t(X_train),R)
      train = cbind(y_train,XRP_train)
      train1 = train[train[,1]==1,]
      train2 = train[train[,1]==2,]
      mu1 = colMeans(train1[,-1])
      mu2 = colMeans(train2[,-1])
      sig1 = crossprod((train1[,-1]-mu1),(train1[,-1]-mu1))/(n1-1)
      sig2 = crossprod((train2[,-1]-mu2),(train2[,-1]-mu2))/(n2-1)
      q_RP_hat = discriminant(X_test%*%R,mu1,mu2,sig1,sig2,n1,n2)
      return(q_RP_hat)
    }
  }
  else if(type == '+1/-1'){
    q_RP_hat1 = foreach(b = 1:B, .combine = 'cbind', .export = c('discriminant','logdet')) %dopar% {
      r = rbinom(d*p,1,0.5)
      r[r==0] = -1
      R = matrix(r/d,nrow=p,ncol=d)
      XRP_train = crossprod(t(X_train),R)
      train = cbind(y_train,XRP_train)
      train1 = train[train[,1]==1,]
      train2 = train[train[,1]==2,]
      mu1 = colMeans(train1[,-1])
      mu2 = colMeans(train2[,-1])
      sig1 = crossprod((train1[,-1]-mu1),(train1[,-1]-mu1))/(n1-1)
      sig2 = crossprod((train2[,-1]-mu2),(train2[,-1]-mu2))/(n2-1)
      q_RP_hat = discriminant(X_test%*%R,mu1,mu2,sig1,sig2,n1,n2)
      return(q_RP_hat)
    }
  }
  else if(type == '3_point'){
    
  q_RP_hat1 = foreach(b = 1:B, .combine = 'cbind',.export = c('discriminant','logdet'), .packages = c('Matrix','doParallel')) %dopar% {

      R <- matrix(sample(c(1,0,-1),(d*p),replace = TRUE, prob = c(1/(2*sqrt(p)), 1 - 1/sqrt(p), 1/(2*sqrt(p)))), nrow = p, ncol = d)
       XRP_train = sapply(1:d, function(i){
          r.loop = R[,i];
          r.index.minus1 = which(r.loop == -1);
          r.index.1 = which(r.loop == 1);
          xr = rowSums(X_train[ ,r.index.1]) - rowSums(X_train[ ,r.index.minus1]) 
          return(xr)
        })

      train = cbind(y_train,XRP_train)
      
      train1 = train[train[,1] == 1, ]
      train2 = train[train[,1] == 2, ]
      mu1 = (colMeans(train1[,-1]))
      mu2 = (colMeans(train2[,-1]))
      sig1 = tcrossprod((t(train1[,-1])-mu1),(t(train1[,-1])-mu1))/(n1-1)
      sig2 = tcrossprod((t(train2[,-1])-mu2),(t(train2[,-1])-mu2))/(n2-1)
      q_RP_hat <- discriminant(X_test%*%R, mu1, mu2, sig1, sig2, n1, n2)
      return(q_RP_hat)
    }
  }else if(type == 'Haar'){
    
    q_RP_hat <- foreach(b = 1:B, .combine = 'cbind',.export = c('discriminant','logdet'), .packages = c('Matrix','doParallel')) %dopar% {
      
      R = matrix(rnorm(d * p), nrow = p, ncol = d)
      R.svd = svd(R)
      R = R.svd$u
      XRP_train <- crossprod(t(X_train),R)
      train <- cbind(y_train,XRP_train)
      train1 <- train[train[,1]==1,]
      train2 <- train[train[,1]==2,]
      mu1 <- colMeans(train1[,-1])
      mu2 <- colMeans(train2[,-1])
      sig1 <- crossprod((train1[,-1]-mu1),(train1[,-1]-mu1))/(n1-1)
      sig2 <- crossprod((train2[,-1]-mu2),(train2[,-1]-mu2))/(n2-1)
      q_RP_hat <- discriminant(X_test%*%R,mu1,mu2,sig1,sig2,n1,n2)
      return(q_RP_hat)
    }
  }
  dis_RP <- rowMeans(q_RP_hat1)
  class <- ifelse(dis_RP>0,1,2)
  return(class)
}





