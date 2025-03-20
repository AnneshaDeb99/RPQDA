logdet <- function(A) {
  L <- eigen(A)$values
  a <- sum(log(abs(L)))
  return(a)
}
a=0

private_BayesRule <- function(X, Omega, Delta, Logratio, u1, Delta_S2) {
  n <- nrow(X)
  p <- ncol(X)
  X <- t(X)
  results <- diag(t(X - matrix(rep(u1, n), nrow = length(u1), ncol = n)) %*% Omega %*% (X - matrix(rep(u1, n), nrow = length(u1), ncol = n))) - 
    2 * as.vector(t(Delta_S2 %*% (X - matrix(rep(u1, n), nrow = length(u1), ncol = n)))) + as.vector(Delta_S2 %*% Delta) - Logratio
  
  
  results <- ifelse((results > 0),1,0)
  return(results)
}



private_bisectTrainDelta <- function(testSet1, testSet2, Omega, Logratio, u1, u2, invS2) {
  Delta0 <- u2 - u1
  th1 <- 0
  th2 <- max(abs(Delta0))
  n <- nrow(testSet1) + nrow(testSet2)
  
  while (abs(th1 - th2) > 1e-3) {
    Delta <- Delta0 * (abs(Delta0) > th1)
    Delta_S2 <- t(Delta) %*% invS2
    result1 <- private_BayesRule(testSet1, Omega, Delta, Logratio, u1, Delta_S2)
    result2 <- private_BayesRule(testSet2, Omega, Delta, Logratio, u1, Delta_S2)
    misrate <- (sum(result1 != 1) + sum(result2 != 0)) / n
    
    Delta <- Delta0 * (abs(Delta0) > th2)
    Delta_S2 <- t(Delta) %*% invS2
    result1 <- private_BayesRule(testSet1, Omega, Delta, Logratio, u1, Delta_S2)
    result2 <- private_BayesRule(testSet2, Omega, Delta, Logratio, u1, Delta_S2)
    
    if ((sum(result1 != 1) + sum(result2 != 0)) / n < misrate) {
      th1 <- (th1 + th2) / 2
      misrate <- (sum(result1 != 1) + sum(result2 != 0)) / n
    } else {
      th2 <- (th1 + th2) / 2
    }
    
    Delta <- Delta0 * (abs(Delta0) > th2)
  }
  
  return(list(Delta = Delta, th2 = th2, misrate = misrate))
}

private_bisectTrainSQDA <- function(X1, X2, nfolds) {
  n1 <- nrow(X1)
  p <- ncol(X1)
  n2 <- nrow(X2)
  
  test1 <- private_produceValidationSets(n1, n2, nfolds)$test1
  test2 <- private_produceValidationSets(n1, n2, nfolds)$test2
  
  u1 <- apply(X1, 2, mean)
  u2 <- apply(X2, 2, mean)
  S1 <- cov(X1)
  S2 <- cov(X2)
  
  Delta0 <- u2 - u1
  Sdiff <- S2 - S1
  Thresholds <- matrix(0, 3, 2)
  
  Thresholds[2, 1] <- quantile(abs(Sdiff), 0.5)
  Thresholds[2, 2] <- max(abs(Sdiff))
  
  Thresholds[3, 1] <- quantile(abs(S2[!diag(nrow(S2))]), 0.5)
  Thresholds[3, 2] <- max(abs(S2[!diag(nrow(S2))]))
  
  Thresholds[1, 1] <- 0
  Thresholds[1, 2] <- max(abs(Delta0))
  
  while (min(abs(Thresholds[3, 2] - Thresholds[3, 1]), abs(Thresholds[2, 2] - Thresholds[2, 1]), abs(Thresholds[2, 2] - Thresholds[2, 1])) > 1e-4) {
    results <- array(0, c(2, 2, 2))
    for (t in 1:nfolds) {
      testSet1 <- X1[test1[[t]], ]
      testSet2 <- X2[test2[[t]], ]
      testSize <- length(test1[[t]]) + length(test2[[t]])
      trainSet1 <- X1[-test1[[t]], ]
      trainSet2 <- X2[-test2[[t]], ]
      
      u1 <- apply(trainSet1, 2, mean)
      u2 <- apply(trainSet2, 2, mean)
      
      S1 <- cov(trainSet1)
      S2 <- cov(trainSet2)
      
      for (i in 1:2) {
        for (j in 1:2) {
          for (k in 1:2) {
            th1 <- Thresholds[1, i]
            th2 <- Thresholds[2, j]
            th3 <- Thresholds[3, k]
            result <- private_thresholding(S1, S2, u1, u2, th1, th2, th3)
            Omega <- result$Omega
            Delta <- result$Delta
            Logratio <- result$Logratio
            Delta_S2 <- result$Delta_S2
            result1 <- private_BayesRule(testSet1, Omega, Delta, Logratio, u1, Delta_S2)
            result2 <- private_BayesRule(testSet2, Omega, Delta, Logratio, u1, Delta_S2)
            results[i, j, k] <- results[i, j, k] + ((sum(result1 != 1) + sum(result2 != 0)) / testSize)
          }
        }
      }
    }
    
    # Find the minimum
    misrate <- 100
    for (i in 1:2) {
      for (j in 1:2) {
        for (k in 1:2) {
          if (results[i, j, k] < misrate) {
            th <- c(i, j, k)
            misrate <- results[i, j, k]
          }
        }
      }
    }
    
    for (i in 1:3) {
      Thresholds[i, 3 - th[i]] <- mean(Thresholds[i, ])
    }
  }
  
  u1 <- apply(X1, 2, mean)
  u2 <- apply(X2, 2, mean)
  S1 <- cov(X1)
  S2 <- cov(X2)
  
  result <- private_thresholding(S1, S2, u1, u2, Thresholds[1, 2], Thresholds[2, 2], Thresholds[3, 2])
  Omega <- result$Omega
  Delta <- result$Delta
  Logratio <- result$Logratio
  Delta_S2 <- result$Delta_S2
  
  return(list(Omega = Omega, Delta = Delta, Logratio = Logratio, u1 = u1, Delta_S2 = Delta_S2))
}



private_produceValidationSets <- function(n1, n2, nfolds) {
  v1 <- floor(n1 / nfolds)
  v2 <- floor(n2 / nfolds)
  
  test1 <- vector("list", nfolds)
  test2 <- vector("list", nfolds)
  
  for (i in 1:(nfolds - 1)) {
    test1[[i]] <- ((1 + (i - 1) * v1):(i * v1))
    test2[[i]] <- ((1 + (i - 1) * v2):(i * v2))
  }
  
  test1[[nfolds]] <- ((1 + (nfolds - 1) * v1):n1)
  test2[[nfolds]] <- ((1 + (nfolds - 1) * v2):n2)
  
  return(list(test1 = test1, test2 = test2))
}


private_quantileTrainSQDA <- function(X1, X2, nfolds, nh1, nh2, nh3) {
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  p <- ncol(X1)
  
  # Produce the cross-validation set
  validation_sets <- private_produceValidationSets(n1, n2, nfolds)
  test1 <- validation_sets$test1
  test2 <- validation_sets$test2
  
  # Determine the threshold values
  u1 <- colMeans(X1)
  u2 <- colMeans(X2)
  S1 <- cov(X1)
  S2 <- cov(X2)
  
  Delta0 <- u2 - u1
  Sdiff <- S2 - S1
  
  dh1 <- 50 / nh1
  qh1 <- seq(50, 100, by = dh1) / 100
  h1 <- quantile(abs(Delta0), qh1)
  
  dh2 <- 50 / nh2
  qh2 <- seq(50, 100, by = dh2) / 100
  h2 <- quantile(abs(Sdiff), qh2)
  
  dh3 <- 50 / nh3
  qh3 <- seq(50, 100, by = dh3) / 100
  h3 <- quantile(abs(S2[!diag(nrow(S2))]), qh3)
  
  results <- array(0, dim = c(length(h1), length(h2), length(h3)))
  
  # Start cross-validation
  for (t in 1:nfolds) {
    testSet1 <- X1[test1[[t]], ]
    testSet2 <- X2[test2[[t]], ]
    testSize <- length(test1[[t]]) + length(test2[[t]])
    
    trainSet1 <- X1[-test1[[t]], ]
    trainSet2 <- X2[-test2[[t]], ]
    
    u1 <- colMeans(trainSet1)
    u2 <- colMeans(trainSet2)
    
    S1 <- cov(trainSet1)
    S2 <- cov(trainSet2)
    
    for (i in 1:length(h1)) {
      for (j in 1:length(h2)) {
        for (k in 1:length(h3)) {
          thresholding_result <- private_thresholding(S1, S2, u1, u2, h1[i], h2[j], h3[k])
          Omega <- thresholding_result$Omega
          Delta <- thresholding_result$Delta
          Logratio <- thresholding_result$Logratio
          Delta_S2 <- thresholding_result$Delta_S2
          
          result1 <- private_BayesRule(testSet1, Omega, Delta, Logratio, u1, Delta_S2)
          result2 <- private_BayesRule(testSet2, Omega, Delta, Logratio, u1, Delta_S2)
          
          results[i, j, k] <- results[i, j, k] + ((sum(result1 != 1) + sum(result2 != 0)) / testSize)
        }
      }
      cat("Thresholds ", h1[i], " for fold ", t, " finished!\n")
    }
  }
  
  # Find the minimum
  misrate <- 100
  for (i in 1:length(h1)) {
    for (j in 1:length(h2)) {
      for (k in 1:length(h3)) {
        
        if (results[i, j, k] < misrate) {
          th <- c(i, j, k)
          misrate <- results[i, j, k]
        }
      }
    }
  }
  
  u1 <- colMeans(X1)
  u2 <- colMeans(X2)
  S1 <- cov(X1)
  S2 <- cov(X2)
  
  final_result <- private_thresholding(S1, S2, u1, u2, h1[th[1]], h2[th[2]], h3[th[3]])
  Omega <- final_result$Omega
  Delta <- final_result$Delta
  Logratio <- final_result$Logratio
  Delta_S2 <- final_result$Delta_S2
  
  return(list(Omega = Omega, Delta = Delta, Logratio = Logratio, u1 = u1, Delta_S2 = Delta_S2))
}




private_quantileTrainSQDA2 <- function(X1, X2, nfolds, nh1, nh2, nh3) {
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  p <- ncol(X1)
  
  # Produce the cross-validation set
  validation_sets <- private_produceValidationSets(n1, n2, nfolds)
  test1 <- validation_sets$test1
  test2 <- validation_sets$test2
  
  # Determine the threshold values
  u1 <- colMeans(X1)
  u2 <- colMeans(X2)
  S1 <- cov(X1)
  S2 <- cov(X2)
  
  Sdiff <- S2 - S1
  
  # dh1 <- 50 / nh1
  # qh1 <- seq(50, 100, by = dh1) / 100
  # h1 <- quantile(abs(Delta0), qh1)
  
  dh2 <- 50 / nh2
  qh2 <- seq(50, 100, by = dh2) / 100
  h2 <- quantile(abs(Sdiff), qh2)
  
  dh3 <- 50 / nh3
  qh3 <- seq(50, 100, by = dh3) / 100
  h3 <- quantile(abs(S2[!diag(nrow(S2))]), qh3)
  
  results <- matrix(0, nrow = length(h2), ncol = length(h3))
  th <- matrix(0, nrow = length(h2), ncol = length(h3))
  
  # Start cross-validation
  for (t in 1:nfolds) {
    testSet1 <- X1[test1[[t]], ]
    testSet2 <- X2[test2[[t]], ]
    testSize <- length(test1[[t]]) + length(test2[[t]])
    
    trainSet1 <- X1[-test1[[t]], ]
    trainSet2 <- X2[-test2[[t]], ]
    
    u1 <- colMeans(trainSet1)
    u2 <- colMeans(trainSet2)
    
    S1 <- cov(trainSet1)
    S2 <- cov(trainSet2)
    
    for (j in 1:length(h2)) {
      for (k in 1:length(h3)) {
        thresholding_result <- private_thresholding(S1, S2, u1, u2, 0, h2[j], h3[k])
        Omega <- thresholding_result$Omega
        Logratio <- thresholding_result$Logratio
        invS2 <- thresholding_result$invS2
        
        bisect_result <- private_bisectTrainDelta(testSet1, testSet2, Omega, Logratio, u1, u2, invS2)
        Delta <- bisect_result$Delta
        th2 <- bisect_result$th2
        misrate <- bisect_result$misrate
        
        results[j, k] <- results[j, k] + misrate
        th[j, k] <- th[j, k] + th2
        
        cat("Thresholds ", h2[j], " for fold ", t, " finished!\n")
      }
    }
  }
  th <- th / nfolds
  
  # Find the minimum
  misrate <- 100
  for (j in 1:length(h2)) {
    for (k in 1:length(h3)) {
      if (results[j, k] < misrate) {
        thresholds <- c(th[j, k], h2[j], h3[k])
        misrate <- results[j, k]
      }
    }
  }
  
  u1 <- colMeans(X1)
  u2 <- colMeans(X2)
  S1 <- cov(X1)
  S2 <- cov(X2)
  
  final_result <- private_thresholding(S1, S2, u1, u2, thresholds[1], thresholds[2], thresholds[3])
  Omega <- final_result$Omega
  Delta <- final_result$Delta
  Logratio <- final_result$Logratio
  Delta_S2 <- final_result$Delta_S2
  
  return(list(Omega = Omega, Delta = Delta, Logratio = Logratio, u1 = u1, Delta_S2 = Delta_S2))
}



private_thresholding <- function(S1, S2, u1, u2, th1, th2, th3) {
  Delta <- u2 - u1
  Sdiff <- S2 - S1
  p <- nrow(Sdiff)
  
  Delta <- Delta * (abs(Delta) > th1)
  Sdiff <- Sdiff * (abs(Sdiff) > th2)
  
  S1[Sdiff == 0] <- S1[Sdiff == 0] + S2[Sdiff == 0] / 2
  S2[Sdiff == 0] <- S1[Sdiff == 0]
  
  diag1 <- diag(S1)
  diag2 <- diag(S2)
  
  S1 <- S1 * (abs(S1) > th3)
  S2 <- S2 * (abs(S2) > th3)
  
  diag(S1) <- diag1
  diag(S2) <- diag2
  
  invS1 <- ginv(S1)
  invS2 <- ginv(S2)
  
  Omega <- invS2 - invS1
  Logratio <- logdet(S1) - logdet(S2)
  Delta_S2 <- t(Delta) %*% invS2
  
  return(list(Omega = Omega, Delta = Delta, Logratio = Logratio, Delta_S2 = Delta_S2, invS2 = invS2))
}



DAQDA <- function(X1, X2, options = list(train = 'bisect',nfolds = 5)) {
  # Set default values
  train <- "bisect"
  nfolds <- 5
  nh1 <- 50
  nh2 <- 25
  nh3 <- 25
  
  # Override defaults with options if provided
  if (!is.null(options$train)) train <- options$train
  if (!is.null(options$nfolds)) nfolds <- options$nfolds
  if (!is.null(options$nDelta)) nh1 <- options$nDelta
  if (!is.null(options$nOmega)) nh2 <- options$nOmega
  if (!is.null(options$nSigma)) nh3 <- options$nSigma
  
  # Training phase
  if (train == "quantile2") {
    results <- private_quantileTrainSQDA2(X1, X2, nfolds, nh1, nh2, nh3)
  } else if (train == "quantile") {
    results <- private_quantileTrainSQDA(X1, X2, nfolds, nh1, nh2, nh3)
  } else if (train == "bisect") {
    results <- private_bisectTrainSQDA(X1, X2, nfolds)
  }
  
  Omega <- results$Omega
  Delta <- results$Delta
  Logratio <- results$Logratio
  u1 <- results$u1
  Delta_S2 <- results$Delta_S2
  
  MODEL <- list(
    Omega = Omega,
    Delta = Delta,
    Logratio = Logratio,
    u1 = u1,
    Delta_S2 = Delta_S2,
    main = which(Delta != 0)
  )
  
  return(MODEL)
}



DAQDAClassify <- function(MODEL, X) {
  Omega <- MODEL$Omega
  Delta <- MODEL$Delta
  Logratio <- MODEL$Logratio
  u1 <- MODEL$u1
  Delta_S2 <- MODEL$Delta_S2
  
  results <- private_BayesRule(X, Omega, Delta, Logratio, u1, Delta_S2)
  return(results)
}

