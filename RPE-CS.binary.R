library(RPEnsemble)

RPBag.binary=function(X_train,y_train,X_test,d=10,B1=500,B2=50){
  Out <- RPParallel(XTrain = X_train, YTrain = (y_train+1), XTest = X_test, d = d,
                    B1 = B1, B2 = B2, base = "QDA", projmethod = "Gaussian", estmethod = "loo",  
                    splitsample = FALSE, k = seq(1, 25, by = 3), clustertype = "Default")
  alphahat <- RPalpha(RP.out = Out, Y = (y_train+1), p1 = 0.5)
  n=length(y_train)
  n.test=length(X_test[,1])
  Class <- RPEnsembleClass(RP.out = Out, n = n,
                           n.test = n.test, p1 = 0.5,alpha = alphahat)
  return(Class)
}
