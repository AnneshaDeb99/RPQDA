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
#### **HDDA**
- **Inputs**:
  - `X_train`: \( n x p \) training data matrix.
  - `y_train`: training classes.
  - `X_test`: \( m x p \) test data matrix.
  - `model`: a parameter of the `hdda` function.
- **Outputs**:
  - `y_hat`: estimated classes.

### 2. `Aoshima-Yata.R`
- **Key Function**: `AoYa`
  #### **AoYa**
- **Inputs**:
  - `X_train`: \( n x p \) training data matrix.
  - `y_train`: training classes.
  - `X_test`: \( m x p \) test data matrix.
- **Outputs**:
  - `y_hat`: estimated classes.

### 3. `IIS-SQDA.R`
- **Key Functions**: `IIS_SQDA`, `predictSQDA`
#### **IIS_SQDA**
- **Inputs**:
  - `X1`: \( n1 x p \) training data matrix for population 1.
  - `X2`: \( n2 x p \) training data matrix for population 2.
  - `options`: parameters for DAQDA.
- **Outputs**:
  - `fit`: the fitted model

#### **predictSQDA**
- **Inputs**:
  - `fit`: the fitted model from `IIS_SQDA` function.
  - `X_test`: \( m \times p \) test data matrix.

- **Outputs**:
  - `y_hat`: estimated classes

### 4. `DA-QDA.R`
- **Key Functions**: `DAQDA`, `DAQDAClassify`

#### **DAQDA**
- **Inputs**:
  - `X1`: \( n1 x p \) training data matrix for population 1.
  - `X2`: \( n2 x p \) training data matrix for population 2.
  - `options`: parameters for DAQDA.
- **Outputs**:
  - `fit`: the fitted model.

#### **DAQDAClassify**
- **Inputs**:
  - `fit`: the fitted model from `DAQDA` function.
  - `X_test`: \( m \times p \) test data matrix.

- **Outputs**:
  - `y_hat`: estimated classes.
  
### 5. `RPE-CS.R`
- **Key Function**: `RPE-CS`
- **Inputs**:
  - `X_train`: \( n \times p \) training data matrix.
  - `y_train`: training classes.
  - `X_test`: \( m \times p \) test data matrix.
  - `d`: reduced dimension.
  - `B1`: number of blocks of matrices.
  - `B2`: number of blocks.
- **Outputs**:
  - `y_hat`: estimated classes.

### 6. `RPE.R`
- **Key Function**: `RPE`
- **Inputs**:
  -`type`: type of random matrix distribution (for example: 'StdNormal' for standard normal)
  - `X_train`: \( n \times p \) training data matrix
  - `y_train`: training classes
  - `X_test`: \( m \times p \) test data matrix
  - `d`: reduced dimension
- **Outputs**:
  - `y_hat`: estimated classes

## Scripts for Model Evaluation
To evaluate the misclassification probabilities for each method and model, the following scripts are used:

### 1. `Model.R`
- **Inputs**:
  - `type`: Scheme/Model type (e.g., `'model1'`, `'model2'`, `'model3'`, `'model4'`)
  - `p`: Data dimension
- **Outputs**:
  - `mu1`, `mu2`: mean paramaters for the Gaussina distribution for class 1 and class 2
  - `sigma1`, `sigma2`: covariance matrices for the Gaussina distribution for class 1 and class 2
  - `Omega1`, `Omega2`: inverse of the covariance matrices for class 1 and class 2
  - `logdet`: logdet(sigma2) - logdet(sigma1)

### 2. `generateData_Bayes.R`
- **Inputs**:
  - `mu1`, `mu2`: mean paramaters for the Gaussina distribution for class 1 and class 2
  - `sigma1`, `sigma2`: covariance matrices for the Gaussina distribution for class 1 and class 2
  - `Omega1`, `Omega2`: inverse of the covariance matrices for class 1 and class 2
  - `logdet`: logdet(sigma2) - logdet(sigma1)
  - `n1`, `n2`: training observations for class 1 and class 2
  -  `m1`, `m2`: test observations for class 1 and class 2
- **Outputs**:
  - `X1`: \( (n1 + m1) x p \) data matrix from class 1
  - `X2`: \( (n2 + m2) x p \) data matrix from class 2
  - `res`: bayes decision rule

### 3. `Experiment.R`
- **Inputs**:
  - `model`: Schemes mentioned in the paper (e.g., `'model1'`, `'model2'`, `'model3'`, `'model4'`)
  - `method`: Methods mentioned in the paper (e.g., `'HDDA'`)
  - `p_all`: vector of dimensions
  - `iter`: number of times data generations to be executed
  - `n1`, `n2`: training observations for class 1 and class 2
  - `m1`, `m2`: test observations for class 1 and class 2
  - `file_name`: file_name where outputs to be stored (e.g., `'result.RData'`)
- **Outputs**:
  - `MP`: list of misclassification proportions of dimension (iter, (n2 + m2), p_all)
  - `Ti`: list of time of execution of dimension (iter, (n2 + m2), p_all)

### 4. `Main_Parallel.R`
This script runs the full evaluation in a parallel computing environment.


