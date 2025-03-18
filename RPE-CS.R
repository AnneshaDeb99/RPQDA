library(RPEnsemble)

class.set = function(j){
  subset = c(which(y_train == set[1,j]), which(y_train == set[2,j]))
  y = as.numeric(as.factor(y_train[subset]))
  X = X_train[subset, ]
  Out = RPParallel(XTrain = X, YTrain = y, XTest = (as.matrix((X_test))), d = d,
                   B1 = B1, B2 = B2, base = "QDA", projmethod = "Gaussian", estmethod = "loo",  
                   splitsample = FALSE, k = seq(1, 25, by = 3), clustertype = "Default")
  #alphahat <- RPalpha(RP.out = Out, Y = y, p1 = 0.5);
  n <- length(y)
  #n.test <- length(X_test[,1])
  n.test <- nrow(as.matrix((X_test)))
  class.subset <- RPEnsembleClass(RP.out = Out, n = n,
                                  n.test = n.test, alpha = RPalpha(Out, Y = y, p1 = sum(y == 1)/length(y)));
  class.subset = ifelse(class.subset == 1, min(set[1,j],set[2,j]), max(set[2,j],set[1,j]))
  return(class.subset)
}

RPBag = function(X_train, y_train, X_test, d = 10, B1 = 500, B2 = 50){
  k = length(unique(y_train))
  c = unique(y_train)
  set = combn(c, 2)
  
  Class <- foreach(j = 1:ncol(set), .combine = c, .export = c('set', 'class.set', 'y_train', 'X_train', 'X_test', 'd', 'B1', 'B2'), .packages = "RPEnsemble")%dopar%{class.set(j)} 
  return(as.numeric(names(which.max(table(Class)))))
}



# function(j){
#     subset = c(which(y_train == set[1,j]), which(y_train == set[2,j]));
#     y = as.numeric(as.factor(y_train[subset]));
#     X = X_train[subset, ];
#     Out = RPParallel(XTrain = X, YTrain = y, XTest = t(as.matrix((X_test))), d = d,
#                      B1 = B1, B2 = B2, base = "QDA", projmethod = "Gaussian", estmethod = "loo",  
#                      splitsample = FALSE, k = seq(1, 25, by = 3), clustertype = "Default");
#     #alphahat <- RPalpha(RP.out = Out, Y = y, p1 = 0.5);
#     n <- length(y);
#     #n.test <- length(X_test[,1])
#     n.test <- nrow(as.matrix(t(X_test)));
#     class.subset <- RPEnsembleClass(RP.out = Out, n = n,
#                                     n.test = n.test, alpha = RPalpha(Out, Y = y, p1 = sum(y == 1)/length(y)));
#     class.subset = ifelse(class.subset == 1, min(set[1,j],set[2,j]), max(set[2,j],set[1,j]))
#     return(class.subset)
