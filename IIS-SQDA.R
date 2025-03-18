# Load required libraries
library(glmnet)
library(clime)
library(glasso)


# Omega_estimate <- function(X){
#   #S <- cov(X)
#   CLIME <- clime(X, standardize = FALSE)
#   lambda <- cv.clime(CLIME)$lambdaopt
#   Omega <- clime(X,lambda = lambda, standardize = FALSE)$Omegalist[[1]]
#   return(Omega)
# }

IIS_SQDA <- function(X1, X2, d = NULL, refit = FALSE) {
  if (is.null(d)) {
    d <- NA
  }
  # Omega1_est <- Omega_estimate(X1)
  # Omega2_est <- Omega_estimate(X2)
  
  # select1 <- IIS(X1, X2, Omega1_est, d)
  # select2 <- IIS(X1, X2, Omega2_est, d)
  
  Omega1 <- glasso(cov(X1),rho=0.1)$wi
  Omega2 <- glasso(cov(X2),rho=0.1)$wi
  
  select1 <- IIS(X1, X2, Omega1, d)
  select2 <- IIS(X1, X2, Omega2, d)
  select <- unique(c(select1, select2))
  
  # Combine the datasets
  X <- rbind(X1, X2)
  Y <- c(rep(1, nrow(X1)), rep(0, nrow(X2)))
  
  # Generate interaction terms
  Xintr2 <- NULL
  for (i in 1:length(select)) {
    for (j in i:length(select)) {
      Xintr2 <- cbind(Xintr2, X[, select[i]] * X[, select[j]])
    }
  }
  X2 <- cbind(X, Xintr2)
  
  # Logistic regression with cross-validation
  
  output1 <- cv.glmnet(X2, Y, family = "binomial", alpha = 0.5)
  support <- sort(which(as.vector(coef(output1, s = "lambda.min"))[-1] != 0))
  X2New <- X2[, support, drop = FALSE]
  
  if (refit) {
    model <- glm(Y ~ ., data = as.data.frame(X2New), family = binomial(link = "logit"))
  } else {
    model <- output1
  }
  
  output <- list()
  output$model <- model
  output$refit <- refit
  output$main <- which(support <= ncol(X))
  output$support <- support
  output$select <- select
  
  support <- setdiff(support, 1:ncol(X))
  inter_matrix <- matrix(0, nrow = ncol(X), ncol = ncol(X))
  current_location <- ncol(X)
  
  for (i in 1:length(select)) {
    for (j in i:length(select)) {
      current_location <- current_location + 1
      if (length(support) == 0) {
        break
      } else if (current_location == support[1]) {
        inter_matrix[select[i], select[j]] <- 1
        inter_matrix[select[j], select[i]] <- 1
        support <- support[-1]
      }
    }
  }
  
  output$inter <- inter_matrix
  return(output)
}

# Example call to the function
# Note: The IIS function needs to be defined in R
# result <- IIS_SQDA(X1, X2, Omega1, Omega2, d, refit)

IIS <- function(X1, X2, Omega, d = NULL) {
  # Compute the transformed variables
  Z1 <- as.matrix(X1) %*% as.matrix(Omega)
  Z2 <- as.matrix(X2) %*% as.matrix(Omega)
  
  # Determine the number of samples
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  
  # Determine d if not provided
  if (is.null(d)) {
    d <- max(floor((n1 + n2) / 10), 10)
  }
  
  # Compute the standard deviations
  sigmaZ1 <- apply(Z1, 2, sd)
  sigmaZ2 <- apply(Z2, 2, sd)
  sigmaZ12 <- apply(rbind(Z1, Z2), 2, sd)
  
  # Compute the test statistic D
  D1 <- log(sigmaZ12^2) - (n1 / (n1 + n2)) * log(sigmaZ1^2) - (n2 / (n1 + n2)) * log(sigmaZ2^2)
  
  # Sort and select the top d indices
  index <- order(D1, decreasing = TRUE)
  output <- index[1:d]
  
  return(output)
}

# Example call to the function
# Note: Ensure X1, X2, and Omega are defined as matrices
# result <- IIS(X1, X2, Omega, d)
predictSQDA <- function(model, newX) {
  library(glmnet)
  
  Xintr2 <- newX
  select <- model$select
  support <- model$support
  
  # Create interaction terms
  for (i in seq_along(select)) {
    for (j in i:length(select)) {
      Xintr2 <- cbind(Xintr2, newX[, select[i]] * newX[, select[j]])
    }
  }
  
  X2new <- Xintr2[, support]
  
  if (model$refit) {
    predictions <- predict(model$model, newdata = as.data.frame(X2new), type = "response")
  } else {
    predictions <- predict(model$model, newx = as.matrix(Xintr2), s = "lambda.min", type = "response")
  }
  
  output <- ifelse((predictions > 0.5),0,1)
  return(output)
}

# Example call to the function
# Note: Ensure `model` is defined and `newX` is a matrix or data frame
# result <- predictSQDA(model, newX)



