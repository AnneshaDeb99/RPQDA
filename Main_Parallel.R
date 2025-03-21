################################################################################
################################################################################

# 1. Model.R will give us population parameters mu1, mu2, sigma1 and sigma2
#
# 2. Given the parameters generateData_Bayes.R will generate data (X1) and (X2)
#    and calculates class assignments using bayes rule (res)
#
# 3. Using the generated data, Experiment.R will calculate MP (Misclassification
#    proportion) and Ti (Time) for each methods (e.g. Bayes, RPE-SN etc.), for
#    each dimensions and iterations.
#
################################################################################
################################################################################

rm(list = ls())

#Provide the path where other files are stored.
path <- readline("Enter the path to the Codes folder: ")
setwd(path)

library(MASS)
library(parallel)
library(doSNOW)
library(doParallel)
  

source('Experiment.R')

cl <- makeCluster(8, type = "SOCK")
registerDoSNOW(cl)

scheme <- c('scheme1','scheme2','scheme3','scheme4')  #Specifying the model, types of model are specified in Model.R file


method <- 'all'

p_all <- c(2^(8:13),10000)   #all the dimensions to be evaluated

p_all = 100

n1 <- 100 #train sample size population 1
n2 <- 100 #train sample size population 2
m1 <- 200 #test sample size population 1
m2 <- 200 #test sample size population 2

iter <- 10

run <- function(i){
  file_name <- paste0(scheme[i], 'of', method, '.RData')
  cat("Starting model", i, "...\n")
  result <- Experiment(scheme = scheme[i],method = method, p_all = p_all,iter = iter,n1 = n1,n2 = n2, m1 = m1, m2 = m2, file_name = file_name)
  MP <- result$MP
  Ti <- result$Ti
  save(MP,Ti, file = file_name)
}


print(date())
foreach(i = 1:length(scheme), .packages = 'doParallel')%dopar%{run(i)}
print(date())



