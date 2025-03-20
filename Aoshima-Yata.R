library(psych)


discrim = function(j, y_train, X_train, x.te){
  n = sum(y_train == j);
  class = which(y_train == j);
  xay.train = as.matrix(X_train[class, ]);
  mu = colMeans(xay.train);
  sig = crossprod(xay.train - mu, xay.train - mu)/n;
  if(ncol(x.te - mu) == 1){
    dis = crossprod(x.te - mu, x.te - mu) - (tr(sig)/n);
  }else{
    dis = tcrossprod(x.te - mu, x.te - mu) - (tr(sig)/n);
  }
  return(dis)
}

#

AoYa=function(X_train, y_train, X_test){
  
  if(is.vector(X_test) == TRUE){
    X_test = t(as.matrix(X_test))
  }
  k <- length(unique(y_train)) #number of classes
  p <- ncol(X_test)
  n.test <- nrow(X_test)
  q.AY.hat = matrix(0, nrow = n.test, ncol = k)
  class = NULL
  #y_train = ifelse(y_train == 0, 1, 2)
  
  for(i in 1:n.test){
    
    x.te = as.matrix(X_test[i, ])
    
    q.AY.hat[i, ] = foreach(j = 1:k, .combine = c, .packages = 'psych', .export = c('discrim'))%dopar%{discrim(j, y_train, X_train, x.te)}
    class[i] = max(which(q.AY.hat[i, ] == (min(q.AY.hat[i, ]))))
  }
  #class = ifelse(class == 1, 0, 1)
  return(class)
}