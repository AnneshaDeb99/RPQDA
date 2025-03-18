################################################################################
################################################################################

# 1. Model.R will give us population parameters mu1, mu2, sigma1 and sigma2
#
# 2. Given the parameters generateData_Bayes.R will generate data (X1) and (X2)
#    and calculates class assignments using bayes rule (res)
#
# 3. Using the generated data, Experiment.R will calculate EMP (Misclassification
#    proportion) and Ti (Time) for each methods (e.g. Bayes, RPE-CS etc.), for
#    each dimensions and iterations.
#
# 4. Change all the neccessary parameters such as for RPE-SN and RPE-TP , d and
#    B in Experiment.R file.
#

################################################################################
################################################################################

#New
rm(list = ls())

#Working
path <- "C:/Users/Annsha deb/Dropbox/Random_Projection/QDA/QDA_Codes/Codes_for_the_paper"

#path <- readline("Enter the path to the Codes folder: ")
setwd(path)

library(MASS)
library(tictoc)
library(parallel)
library(doSNOW)
library(doParallel)
  

source('Experiment.R')

cl <- makeCluster(8, type = "SOCK")
registerDoSNOW(cl)

model <- c('model1','model2','model4','model5')  #Specifying the model, types of model are specified in Model.R file

#Specify the method : 'all' to evaluate all the methods considered here, or
#                     type the name 'Bayes', 'HDClassif', 'AoYa', 'RPBag', 'sQDA', 'IIS_SQDA', 'RPE_SN', 'RPE_TP'.

method <- 'RPE-TP'

p_all <- c(2^(8:13),10000)   #all the dimensions to be evaluated


n1 <- 100 #train sample size population 1
n2 <- 100 #train sample size population 2
m1 <- 200 #test sample size population 1
m2 <- 200 #test sample size population 2

iter <- 20

run <- function(i){
  file_name <- paste0(model[i], 'of', method, '.RData')
  #tic()
  cat("Starting model", i, "...\n")
  result <- Experiment(model = model[i],method = method, p_all = p_all,iter = iter,n1 = n1,n2 = n2, m1 = m1, m2 = m2, file_name = file_name)
  #toc()
  MP <- result$MP
  Ti <- result$Ti
  save(MP,Ti, file = file_name)
}


print(date())
foreach(i = 1:length(model), .packages = 'doParallel')%dopar%{run(i)}
print(date())



