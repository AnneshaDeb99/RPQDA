library(HDclassif)
HDDA = function(X_train,y_train,X_test, model = "ALL"){
  model = hdda(X_train,y_train,model = model)
  y_hat = predict(model,X_test)$class
  return(y_hat)
}
