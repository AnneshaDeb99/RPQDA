# RPQDA

The following description of the functions refers to the paper titled *"Ultrahigh-dimensional Quadratic Discriminant Analysis Using Random Projections"*.

## Competing Methods
The competing methods include:
1. **HDDA**
2. **IIS-SQDA**
3. **DA-QDA**
4. **RPE-CS** (Cannings-Samworth)
5. **AoYa** (Aoshima-Yata)

Our proposed method: **RPE**

## Function Descriptions

### 1. `HDDA.R` 
- **Key Function**: `HDDA`
- **Inputs**:
  - `X_train`: \( n x p \) training data matrix
  - `y_train`: training classes
  - `X_test`: \( m x p \) test data matrix
  - `model`: a parameter of the `hdda` function
- **Outputs**:
  - `y_hat`: estimated classes

### 2. `Aoshima-Yata.R`
- **Key Function**: `AoYa`
- **Inputs**:
  - `X_train`: \( n x p \) training data matrix
  - `y_train`: training classes
  - `X_test`: \( m x p \) test data matrix
- **Outputs**:
  - `y_hat`: estimated classes

### 3. `IIS-SQDA.R`
- **Key Functions**: `IIS_SQDA`, `predictSQDA`
- **IIS_SQDA**
- **Inputs**:
  - `X1`: \( n1 x p \) training data matrix for population 1
  - `X2`: \( n2 x p \) training data matrix for population 2
  - `options`: parameters for DAQDA
- **Outputs**:
  - `fit`: the fitted model

- **predictSQDA**:
- **Inputs**:
  - `fit`: the fitted model from `IIS_SQDA` function
  - `X_test`: \( m \times p \) test data matrix

- **Outputs**:
  - `y_hat`: estimated classes

### 4. `DA-QDA.R`
- **Key Functions**: `DAQDA`, `DAQDAClassify`
- **DAQDA**:
- **Inputs**:
  - `X1`: \( n1 x p \) training data matrix for population 1
  - `X2`: \( n2 x p \) training data matrix for population 2
  - `options`: parameters for DAQDA
- **Outputs**:
  - `fit`: the fitted model

- **DAQDAClassify**:
- **Inputs**:
  - `fit`: the fitted model from `DAQDA` function
  - `X_test`: \( m \times p \) test data matrix

- **Outputs**:
  - `y_hat`: estimated classes
  
### 5. `RPE-CS.R`
- **Key Function**: `RPE-CS`
- **Inputs**:
  - `X_train`: \( n \times p \) training data matrix
  - `y_train`: training classes
  - `X_test`: \( m \times p \) test data matrix
  - `d`: reduced dimension
  - `B1`: number of blocks of matrices
  - `B2`: number of blocks
- **Outputs**:
  - `y_hat`: estimated classes

### 6. `RPE.R`
- **Key Function**: `RPE`
- **Inputs**:
  - `X_train`: \( n \times p \) training data matrix
  - `y_train`: training classes
  - `X_test`: \( m \times p \) test data matrix
- **Outputs**:
  - `y_hat`: estimated classes

## Scripts for Model Evaluation
To evaluate the misclassification probabilities for each method, the following scripts are used:

### 1. `Model.R`
- **Inputs**:
  - `type`: Model type (e.g., `'model1'`, `'model2'`, `'model3'`, `'model4'`)
  - `p`: Data dimension
- **Outputs**:
  - `mu1`, `mu2`
  - `sigma1`, `sigma2`
  - `Omega1`, `Omega2`
  - `logdet`

### 2. `generateData_Bayes.R`
- **Inputs**:
  - `mu1`, `mu2`
  - `sigma1`, `sigma2`
  - `Omega1`, `Omega2`
  - `logdet`
  - `n1`, `n2`, `m1`, `m2`
- **Outputs**:
  - `X1`
  - `X2`
  - `res`

### 3. `Experiment.R`
- **Inputs**:
  - `mu1`, `mu2`
  - `sigma1`, `sigma2`
  - `Omega1`, `Omega2`
  - `logdet`
  - `n1`, `n2`, `m1`, `m2`
- **Outputs**:
  - `X1`
  - `X2`
  - `res`

### 4. `Main_Parallel.R`
This script runs the full evaluation in a parallel computing environment.

The following description of the functions are in reference to the paper titled "Ultrahigh-dimensional Quadratic Discriminant Analysis Using Random Projections".

The competeting methods are - 1. HDDA, 2. IIS-SQDA, 3. DA-QDA, 4. RPE-CS (Cannings-Samworth), 5. AoYa (Aoshima-Yata).
Our method - RPE

1. HDDA.R : This script is based on the hdda function from the library HDClassif.

   Key function: HDDA

      Inputs: (i) X_train : n x p train data matrix, (ii) y_train: training classes, (iii) X_test: m x p test data matrix, (iv) model: a parameter of hdda function

      Outputs: (i) y_hat: estimated classes 

3. Aoshima-Yata.R : 
   Key function: AoYa
       Inputs: (i) X_train : n x p train data matrix, (ii) y_train: training classes, (iii) X_test: m x p test data matrix.
       Outputs: (i) y_hat: estimated classes.

4. IIS-SQDA.R
   Key function:
      Inputs: (i) X_train : n x p train data matrix, (ii) y_train: training classes, (iii) X_test: m x p test data matrix.
      Outputs: (i) y_hat: estimated classes.

5. DA-QDA.R
   Key function:
      Inputs: (i) X_train : n x p train data matrix, (ii) y_train: training classes, (iii) X_test: m x p test data matrix.
      Outputs: (i) y_hat: estimated classes.

6. RPE-CS.R
   Key function:
      Inputs: (i) X_train : n x p train data matrix, (ii) y_train: training classes, (iii) X_test: m x p test data matrix.
      Outputs: (i) y_hat: estimated classes.

7. RPE.R
   Key function:
      Inputs: (i) X_train : n x p train data matrix, (ii) y_train: training classes, (iii) X_test: m x p test data matrix.
      Outputs: (i) y_hat: estimated classes.

We have the scripts Model.R, generateData_Bayes.R, Experiment.R and Main_Parallel.R to evaluate the misclassification probabilities for each of the methods corresponding to the models mentioned in the paper.

   1. Model.R
      Inputs: (i) type: Model types (for example: 'model1', 'model2', 'model3' or 'model4'), (ii) p: data dimension
      Outputs: (i) mu1, mu2: , (ii) sigma1, sigma2: , (iii) Omega1, Omega2: , (iv) logdet

   2. generateData_Bayes.R
      Inputs: (i) mu1, mu2: , (ii) sigma1, sigma2: , (iii) Omega1, Omega2: , (iv) logdet: , (v) n1, n2, m1, m2
      Outputs: (i) X1:, (ii) X2: , (iii) res

   3. Experiment.R
      Inputs: (i) mu1, mu2: , (ii) sigma1, sigma2: , (iii) Omega1, Omega2: , (iv) logdet: , (v) n1, n2, m1, m2
      Outputs: (i) X1:, (ii) X2: , (iii) res

   4. Main_Parallel.R 
