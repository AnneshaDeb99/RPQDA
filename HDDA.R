library(HDclassif)
HDClass=function(X_train,y_train,X_test, model = "ALL"){
  model=hdda(X_train,y_train,model = model)
  res=predict(model,X_test)$class
  return(res)
}
