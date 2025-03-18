#AoYa for two classes


library(psych)
AoYa.binary=function(X_train,y_train,X_test){
  
  ##Data
  train=cbind(y_train,X_train)
  n1=length(train[train[,1]==0,1])
  n2=length(train[,1])-n1
  train1=train[train[,1]==0,]
  train2=train[train[,1]==1,]
  
  ##Parameter Estimation
  mu1=colMeans(train1[,-1])
  mu2=colMeans(train2[,-1])
  sig1=t(train1[,-1]-mu1)%*%(train1[,-1]-mu1)/(n1-1)
  sig2=t(train2[,-1]-mu2)%*%(train2[,-1]-mu2)/(n2-1)
  
  
  ##Discriminant
  dis=(t(t(X_test)-((mu1+mu2)/2))%*%(mu2-mu1))-(tr(sig1)/(2*n1))+(tr(sig2)/(2*n2))
  m=length(dis)
  class=c(0)
  for(i in 1:m)
  {
    if(dis[i]<0)
    {
      class[i]=0
    }
    else
    {
      class[i]=1
    }
  }
  return(class)
}
