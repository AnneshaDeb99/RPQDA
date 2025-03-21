################################################################################
################################## Scheme.R #####################################
################################################################################

# Inputs :  Scheme type (type) and the dimension (p) 
# Outputs:  Means of the two populations: mu1, mu2
#           Covariance matrices: sigma1, sigma2
#           Inverse of the Cov matrices: Omega1, Omega2
#           logdet = logdet(sigma1) - logdet(sigma2)
# Note: The last two are needed in order to calculate the bayes error.

################################################################################
################################################################################
################################################################################


Scheme <- function(type, p){
  library(MASS)
  library(Matrix)
  library(pracma)
  
  delta <- rep(0, p)
  if(p <= 1000){
  delta[c(1, 2)] <- c(0.6, 0.8)
  }else{
    delta[1:(floor(0.002*p))] <- c(rep(c(0.6,0.8),floor(0.001*p)))
  }
  mu2 <- rep(0, p)
  
  #Cannings-Samworth Paper_Scheme1 
  if (type == 'scheme1'){
    
    if(p <= 1000){
      mu1 <- c(rep(3, 3),rep(0, p-3))
    }else{
      mu1 <- c(rep(3, floor(0.003*p)),rep(0, (p-floor(0.003*p))))
    }
    
################################################################################
    
    A <- matrix(rnorm(4 * p), nrow = p, ncol = 4) #(as there are 4 unique eigen values)
    
    # Perform QR decomposition
    ortho <- qr.Q(qr(A))
    
    sigma1 <- ortho%*%diag(c((1+(p-4)*0.5)-0.5,3-0.5,1.5-0.5,1.5-0.5))%*%t(ortho) + 0.5*diag(p)
    Omega1 <- ortho%*%diag(c(1/(1+(p-4)*0.5)-1/0.5,1/3-1/0.5,1/1.5-1/0.5,1/1.5-1/0.5))%*%t(ortho) + (1/0.5)*diag(p)

################################################################################
    
    sigma2 <- ortho[,1:2]%*%diag(c((1+(p-4)*0.5)-0.5,2-0.5))%*%t(ortho[,1:2]) + 0.5*diag(p)
    Omega2 <- ortho[,1:2]%*%diag(c(1/(1+(p-4)*0.5)-1/0.5,1/2-1/0.5))%*%t(ortho[,1:2]) + (1/0.5)*diag(p)
    
################################################################################
    
    logdet1 <- log((2+2*0.5)) + 2*log(2-0.5) + (p-4)*log(1-0.5) + log(1+0.5*(p-4))
    logdet2 <- log(2) + 2*log(0.5) + (p-4)*log(1-0.5) + log(1+0.5*(p-4))
    
    logdet <- 2*logdet2 - 2*logdet1
    
  }
  # Ayoshima-Yata Paper_Scheme2 
  else if (type == 'scheme2'){
  
    #mu1
    mu1 <- rep(0,p)
 
################################################################################    
       
    #mu2
    one <- rep(1,floor(p^(3/5)/2))
    neg_one <- rep(-1,floor(p^(3/5)/2))
    zero <- rep(0,p-2*floor(p^(3/5)/2))
    mu2 <- c(zero,one,neg_one)
    
###############################################################################
        
    #sigma1
    p11 <- floor(p^(2/3))
    p12 <- floor(p^(1/3))
    p13 <- p-p11-p12
    
    rho <- 0.7
    
    sig11 <- (diag(p11)+as.matrix(rep(1,p11))%*%t(rep(1,p11)))/2
    sig12 <- (diag(p12)+as.matrix(rep(1,p12))%*%t(rep(1,p12)))/2
    
    sig11_inv <- (1/0.5)*diag(p11) - matrix(rep((1/(0.5*(p11+1))),(p11*p11)),ncol=p11,nrow=p11)
    det11 <- ((0.5)^(p11-1))*(1+0.5*(p11-1))
    
    sig12_inv <- (1/0.5)*diag(p12) - matrix(rep((1/(0.5*(p12+1))),(p12*p12)),ncol=p12,nrow=p12)
    det12 <- ((0.5)^(p12-1))*(1+0.5*(p12-1))
    
    c1 <- 1

    B1 <- sqrt(0.5+(1:p13/(p13+1)))
    B1_inv <- (1/sqrt(0.5+(1:p13/(p13+1))))
    
    omega1 <- matrix(0,p13,p13)
    for(i in 1:p13)
    {
      for(j in 1:p13)
      {
        #omega1[i,j] <- rho^((abs(i-j))^(1/3))
        omega1[i,j] <- rho^(abs(i-j))
      }
    }
 
    library(matrixcalc)
    sig13 <- c1*hadamard.prod(B1%*%t(B1),omega1)
    
    sig13_inv <- matrix(0,p13,p13)
    sig13_inv[1, 1] <- 1
    for (i in 2:p13) {
      sig13_inv[i, i] <- 1 + rho^2
      sig13_inv[i-1, i] <- - rho
      sig13_inv[i, i-1] <- - rho
    }
    sig13_inv[p13, p13] <- 1
    sig13_inv <- sig13_inv / (1 - rho^2)
    
    sig13_inv <- (1/c1)*hadamard.prod(B1_inv%*%t(B1_inv),sig13_inv)
    sigma1 <- bdiag(sig11,sig12,sig13)
    Omega1 <- bdiag(sig11_inv,sig12_inv,sig13_inv)
    
#############################################################################    
    
    #sigma2
    p21 <- floor(p^(1/2))
    p22 <- floor(p^(1/2))
    p23 <- p-p21-p22
    
    
    rho <- 0.7
    
    sig21 <- (diag(p21)+as.matrix(rep(1,p21))%*%t(rep(1,p21)))/2
    sig22 <- (diag(p22)+as.matrix(rep(1,p22))%*%t(rep(1,p22)))/2
    
    sig21_inv <- (1/0.5)*diag(p21) - matrix(rep((1/(0.5*(p21+1))),(p21*p21)),ncol=p21,nrow=p21)
    det21 <- ((0.5)^(p21-1))*(1+0.5*(p21-1))
    
    sig22_inv <- (1/0.5)*diag(p22) - matrix(rep((1/(0.5*(p22+1))),(p22*p22)),ncol=p22,nrow=p22)
    det22 <- ((0.5)^(p22-1))*(1+0.5*(p22-1))
    
    c2 <- 1.3
    #B2 <- sqrt(0.5+(1:p23/(p23+1)))*diag(p23)
    B2 <- sqrt(0.5+(1:p23/(p23+1)))
    
    #B2_inv <- (1/sqrt(0.5+(1:p23/(p23+1))))*diag(p23)
    B2_inv <- (1/sqrt(0.5+(1:p23/(p23+1))))
    
    omega2 <- matrix(0,p23,p23)
    for(i in 1:p23)
    {
      for(j in 1:p23)
      {
        omega2[i,j] <- rho^(abs(i-j))
      }
    }
    sig23 <- c2*hadamard.prod(B2%*%t(B2),omega2)
    
    sig23_inv <- matrix(0,p23,p23)
    sig23_inv[1, 1] <- 1
    for (i in 2:p23) {
      sig23_inv[i, i] <- 1 + rho^2
      sig23_inv[i-1, i] <- - rho
      sig23_inv[i, i-1] <- - rho
    }
    sig23_inv[p23, p23] <- 1
    sig23_inv <- sig23_inv / (1 - rho^2)
    sig23_inv <- (1/c2)*hadamard.prod(B2_inv%*%t(B2_inv),sig23_inv)
    sigma2 <- bdiag(sig21,sig22,sig23)
  
    Omega2 <- bdiag(sig21_inv,sig22_inv,sig23_inv)
    
################################################################################
    
    logdet <- (p21+p22-p11-p12)*log(0.5) + log((1+0.5*(p21-1))) + log((1+0.5*(p22-1))) - log((1+0.5*(p11-1))) - log((1+0.5*(p12-1))) + (p23-p13)*log(1-rho^2) + p23*log(c2) - p13*log(c1) + 2*sum(log(B2)) - 2*sum(log(B1))
    
    }
   
  else if (type == 'scheme3') {
    
    mu1 <- mu2
    
################################################################################    
    alpha1 <- 0.6
    alpha2 <- 0.7
    beta1 <- 0.4
    beta2 <- 0.3
    rho <- 0.9
    
    p11 <- floor(p^alpha1)
    p21 <- floor(p^alpha2)
    t1 <- floor(p^beta1)
    t2 <- floor(p^beta2)
    p12 <- p - t1*p11
    p22 <- p - t2*p21
    
    one1 <- rep(1,p11)
    one2 <- rep(1,p21)
    equi_mat1 <- c(replicate(t1, (1-rho)*diag(p11) + rho*one1%*%t(one1), simplify = FALSE ), list(diag(p12)))
    equi_mat2 <- c(replicate(t2, (1-rho)*diag(p21) + rho*one2%*%t(one2), simplify = FALSE ), list(diag(p22)))
    
    sigma1 <- bdiag(equi_mat1)
    sigma2 <- bdiag(equi_mat2)
    
    inv_equi1 <- (1/(1-rho))*(diag(p11) - (rho/(1+(p11-1)*rho))*(one1%*%t(one1)))
    inv_equi2 <- (1/(1-rho))*(diag(p21) - (rho/(1+(p21-1)*rho))*one2%*%t(one2))
    
    Omega1 <- bdiag(c(replicate(t1, inv_equi1, simplify = FALSE), list(diag(p12))))
    Omega2 <- bdiag(c(replicate(t2, inv_equi2, simplify = FALSE), list(diag(p22))))
    
################################################################################
    
    logdet <- (p21*t2-t2-p11*t1+t1)*log(1-rho) + t2*log((1+(p21-1)*rho)) - t1*log((1+(p11-1)*rho))
    
  }

  else if (type == 'scheme4') {
    
    mu1 <- mu2
    
################################################################################

      rho <- 0.9
      Omega1 <- diag(p)
      for (i in 1:p) {
        for (j in 1:p) {
          Omega1[i, j] <- rho^abs(i - j)
        }
      }

      c <- 1.5
      Omega2 <- c*Omega1

      sigma1 <- matrix(0,p,p)
      sigma1[1, 1] <- 1
      for (i in 2:p) {
        sigma1[i, i] <- 1 + rho^2
        sigma1[i-1, i] <- - rho
        sigma1[i, i-1] <- - rho
      }
      sigma1[p, p] <- 1
      sigma1 <- sigma1 / (1 - rho^2)

 
    
      sigma2 <- (1/c)*sigma1
  
  ##############################################################################

      logdet <- -p*log(c)

    }
  
    return(list(mu1 = mu1, mu2 = mu2, sigma1 = as.matrix(sigma1), sigma2 = as.matrix(sigma2), Omega1 = as.matrix(Omega1), Omega2 = as.matrix(Omega2), logdet = logdet ))
}
