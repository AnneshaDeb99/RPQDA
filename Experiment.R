################################################################################
############################ Experiment.R ######################################
################################################################################

# Inputs : 
#         1. scheme: scheme type (In Model.R we are generating several schemes)
#         2. method: vector of methods such as Bayes, RPE-SN etc.
#         3. p_all: vector of all the dimensions
#         4. iter: for each dimension and scheme, number of times we want to execute 
#            and calculate the MP (misclassification proportion) and Ti (Time 
#            for execution).
#         5. n1/m1, n2/m2: train/test sample sizes for class 1 and class 2

# Outputs:
#         1. MP: denoting misclassification proportion, it is an array of 
#                dimension 3, 
#                dim 1 represents methods,
#                dim 2 represents iteration (iter)
#                dim 3 represents dimension (p)
#         2. Ti: denoting time for calculating MP, all the other specifications 
#                are the same as MP. 

################################################################################
################################################################################
################################################################################

Experiment <- function(scheme, method, p_all, iter, n1 = 100, n2 = 100, m1 = 200, m2 = 200, file_name) {
  source('Model.R')
  
  MP <- array(0, dim = c(3, iter, length(p_all)))   #Misclassification Proportions
  Ti <- array(0, dim = c(3, iter, length(p_all)))   #Time
  dimnames(MP)[[1]] <- c("Bayes", "RPE-SN", "RPE-TP")
  dimnames(MP)[[3]] <- as.factor(p_all)
  dimnames(Ti)[[1]] <- c("Bayes", "RPE-SN", "RPE-TP")
  dimnames(Ti)[[3]] <- as.factor(p_all)
  for (k in 1:length(p_all)) {
    p <- p_all[k]
    # Fixing the parameters
    set.seed(123)
    param <- Model(type = scheme, p)
  
  for (j in 1:iter) {
    keep <- c("Model", "k","p","param","j","Bayes","Experiment","generateData_Bayes","i","m2","m1","n1","n2","scheme", "path", "Ti", "file_name", "iter","method","MP","p_all","RPEnsemble")
    rm(list = setdiff(ls(), keep))
    source('RPE.R')
    source('generateData_Bayes.R')
    
    set.seed(1000+j)

      dataX <- generateData_Bayes(param$mu1, param$mu2, param$sigma1, param$sigma2, param$Omega1, param$Omega2, param$logdet, n1, n2, m1, m2)
      bayes <- dataX$res
      cat("Starting iteration =", j, ", p =", p, "...\n")

      
      X_train <- rbind(dataX$X1[1:n1,], dataX$X2[1:n2,])
      y_train <- c(rep(1,n1),rep(2,n2))
      X_test <- rbind(dataX$X1[(n1+1):(n1+m1),], dataX$X2[(n2+1):(n2+m2),])
      y_test <- c(rep(1,m1),rep(2,m2))
      
      print(date())
      
      #Bayes
  
      if(method %in% c('Bayes','all'))
      {
        MP[1,j,k] <- mean(y_test!=bayes)
      }
      
      #Our method: RPE-SN or RPE-TP
      if (method %in% c('RPE-SN', 'all')) {
        start_rpe <- Sys.time()
        output <- RPE('StdNormal',X_train,y_train,X_test, d = 10, B = 500)
        MP[2,j,k] <- mean(y_test != (output))
        end_rpe <- Sys.time()
        Ti[2,j,k] <- as.numeric(difftime(end_rpe, start_rpe, units = "secs"))
      }
      
      
      if (method %in% c('RPE-TP', 'all')) {
        start_rpe <- Sys.time()
        output <- RPE('3_point',X_train,y_train,X_test, d = 10, B = 500)
        MP[3,j,k] <- mean(y_test!=output)
        end_rpe <- Sys.time()
        Ti[3,j,k] <- as.numeric(difftime(end_rpe, start_rpe, units = "secs"))
      }

      
      cat("Completed.\n")
      print(date())
      save(MP,Ti, file = file_name)
    }
  }
  
  if(method == 'all'){
    return(list(MP = MP, Ti = Ti))
  }else{
    return(list(MP = MP[method,,], Ti = Ti[method,,])) 
  }
}
  
