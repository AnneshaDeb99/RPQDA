################################################################################
############################ Experiment.R ######################################
################################################################################

# Inputs : 
#         1. model: model type (In Model.R we are generating several models)
#         2. method: vector of methods such as Bayes, RPBag etc.
#         3. p_all: vector of all the dimensions of data
#         4. iter: for each dimens and model, number of times we want to execute 
#            and calculate the MP (misclassification proportion) and Ti (Time 
#            for execution).
#         5. n1, m1, n2, m2: train/test sample sizes for each populations

# Outputs:
#         1. MP: denoting misclassification proportion, it is an array of 
#                dimension 3, 
#                dim 1 represents methods, there are 8 methods including Bayes
#                dim 2 represents iteration (iter)
#                dim 3 represents dimension (p)
#         2. Ti: denoting time for calculating MP, all the other specifications 
#                are the same as MP. 

################################################################################
################################################################################
#New

Experiment <- function(model, method, p_all, iter, n1=100, n2=100, m1=200, m2=200, file_name) {
  source('Model.R')
  
  #Working
 
  
  MP <- array(0, dim = c(8, iter, length(p_all)))   #Misclassification Proportions
  Ti <- array(0, dim = c(8, iter, length(p_all)))   #Time
  dimnames(MP)[[1]] <- c("Bayes","HDDA", "AoYa", "RPE-CS","DA-QDA","IIS","RPE-SN","RPE-TP")
  dimnames(MP)[[3]] <- as.factor(p_all)
  dimnames(Ti)[[1]] <- c("Bayes","HDDA", "AoYa", "RPE-CS","DA-QDA","IIS","RPE-SN","RPE-TP")
  dimnames(Ti)[[3]] <- as.factor(p_all)
  for (k in 1:length(p_all)) {
    p <- p_all[k]
    # Fixing the parameters
    set.seed(123)
    param <- Model(type = model, p)
  
  for (j in 1:iter) {
    keep <- c("Model", "k","p","param","j","Bayes","Experiment","generateData_Bayes","i","IIS_SQDA","m2","m1","n1","n2","model","path","RPBag","sQDA","Ti","AoYa","file_name","HDDA","iter","method","MP","p_all","RPEnsemble")
    rm(list = setdiff(ls(), keep))
    source('HDDA.R')
    source('IIS-SQDA.R')
    source('RPE-Our method.R')
    source('Aoshima-Yata.R')
    #source('AoYa.R')
    source('DA-QDA.R')
    source('generateData_Bayes.R')
    #source('AoYa.binary.R')
    source('RPE-CS.R')
    source('RPE-CS.multiclass.R')
    
    #source('Bayes.R')
    
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

      
      # HDDA
      if (method %in% c('HDDA', 'all')) {
        start_hd <- Sys.time()
        output <- HDDA(X_train, y_train, X_test)
        MP[2,j,k] <- mean(y_test!=output)
        end_hd <- Sys.time()
        Ti[2,j,k] <- as.numeric(difftime(end_hd, start_hd, units = "secs"))
      }
   
      
      # Aoshima-Yata
      if (method %in% c('AoYa', 'all')) {
        start_ay <- Sys.time()
        output <- AoYa(X_train,y_train,X_test)
        MP[3,j,k] <- mean(y_test!=output)
        end_ay <- Sys.time()
        Ti[3,j,k] <- as.numeric(difftime(end_ay, start_ay, units = "secs"))
      }
      
      
      #RPE-CS
      if (method %in% c('RPE-CS', 'all')) {
        start_bag <- Sys.time()
        output <- RPE.CS(X_train,y_train,X_test, d=5, B1 = 500, B2 = 50)
        MP[4,j,k] <- mean(y_test!= output)
        end_bag <- Sys.time()
        Ti[4,j,k] <- as.numeric(difftime(end_bag, start_bag, units = "secs"))
      }
      #print(date())
      
      #DA-QDA
  
      if (method %in% c('DA-QDA','all')) {
        start_sq <- Sys.time()
        #options <- list(train = 'bisect',nfolds = 5)
        fit <- DAQDA(X_train[1:n1,], X_train[(n1+1):(n1+n2),])
        output <- DAQDAClassify(fit, X_test)
        output <- ifelse(output == 1, 1, 2)
        MP[5,j,k] <- mean(y_test!=(1-output))
        end_sq <- Sys.time()
        Ti[5,j,k] <- as.numeric(difftime(end_sq, start_sq, units = "secs"))
      }
     
      
      
      #IIS-SQDA
      if (method %in% c('IIS_SQDA', 'all')) {
        start_iis <- Sys.time()
        fit <- IIS_SQDA(X_train[1:n1,], X_train[(n1+1):(n1+n2),], 30)
        output <- as.vector(predictSQDA(fit, X_test)) + 1
        MP[6,j,k] <- mean(y_test!=output)
        end_iis <- Sys.time()
        Ti[6,j,k] <- as.numeric(difftime(end_iis, start_iis, units = "secs"))
      }


      
      #Our method: RPE-SN or RPE-TP
      if (method %in% c('RPE-SN', 'all')) {
        start_rpe <- Sys.time()
        output <- RPE('StdNormal',X_train,y_train,X_test, d = 10, B = 500)
        MP[7,j,k] <- mean(y_test!=(output))
        end_rpe <- Sys.time()
        Ti[7,j,k] <- as.numeric(difftime(end_rpe, start_rpe, units = "secs"))
      }
      
      
      if (method %in% c('RPE-TP', 'all')) {
        start_rpe <- Sys.time()
        output <- RPE('3_point',X_train,y_train,X_test, d = 10, B = 500)
        MP[8,j,k] <- mean(y_test!=output)
        end_rpe <- Sys.time()
        Ti[8,j,k] <- as.numeric(difftime(end_rpe, start_rpe, units = "secs"))
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
  
