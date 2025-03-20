# RPQDA
The following description of the functions are in reference to the paper titled "Ultrahigh-dimensional Quadratic Discriminant Analysis Using Random Projections".

The competeting methods are - 1. HDDA, 2. IIS-SQDA, 3. DA-QDA, 4. RPE-CS (Cannings-Samworth), 5. AoYa (Aoshima-Yata).
Our method - RPE

1. HDDA.R :

   Inputs:

       (i) X_train
       (ii) y_train
      (iii) X_test
      (iv) model
   Outputs:
   (i) y_hat

   2. Model.R
      Inputs:
      (i) type: Model types (for example: 'model1', 'model2', 'model3' or 'model4')
      (ii) p: data dimension

    Outputs:
    (i) mu1, mu2
     (ii) sigma1, sigma2
     (iii) Omega1, Omega2
     (iv) logdet


   4. generateData_Bayes.R
  
      Inputs:

      (i) mu1, mu2
     (ii) sigma1, sigma2
     (iii) Omega1, Omega2
     (iv) logdet
      (v) n1, n2, m1, m2

      Outputs:

      (i) X1
      (ii) X2
      (iii) res
