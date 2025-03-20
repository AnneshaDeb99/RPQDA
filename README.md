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
  - `X_train`: \( n \times p \) training data matrix.
  - `y_train`: Training classes.
  - `X_test`: \( m \times p \) test data matrix.
  - `model`: A parameter of the `hdda` function.
- **Outputs**:
  - `y_hat`: Estimated classes.

### 2. `Aoshima-Yata.R`
- **Key Function**: `AoYa`
#### **AoYa**
- **Inputs**:
  - `X_train`: \( n \times p \) training data matrix.
  - `y_train`: Training classes.
  - `X_test`: \( m \times p \) test data matrix.
- **Outputs**:
  - `y_hat`: Estimated classes.

### 3. `IIS-SQDA.R`
- **Key Functions**: `IIS_SQDA`, `predictSQDA`
#### **IIS_SQDA**
- **Inputs**:
  - `X1`: \( n1 \times p \) training data matrix for population 1.
  - `X2`: \( n2 \times p \) training data matrix for population 2.
  - `options`: Parameters for DAQDA.
- **Outputs**:
  - `fit`: The fitted model.

#### **predictSQDA**
- **Inputs**:
  - `fit`: The fitted model from `IIS_SQDA` function.
  - `X_test`: \( m \times p \) test data matrix.
- **Outputs**:
  - `y_hat`: Estimated classes.

### 4. `DA-QDA.R`
- **Key Functions**: `DAQDA`, `DAQDAClassify`

#### **DAQDA**
- **Inputs**:
  - `X1`: \( n1 \times p \) training data matrix for population 1.
  - `X2`: \( n2 \times p \) training data matrix for population 2.
  - `options`: Parameters for DAQDA.
- **Outputs**:
  - `fit`: The fitted model.

#### **DAQDAClassify**
- **Inputs**:
  - `fit`: The fitted model from `DAQDA` function.
  - `X_test`: \( m \times p \) test data matrix.
- **Outputs**:
  - `y_hat`: Estimated classes.

### 5. `RPE-CS.R`
- **Key Function**: `RPE.CS`
- **Inputs**:
  - `X_train`: \( n \times p \) training data matrix.
  - `y_train`: Training classes.
  - `X_test`: \( m \times p \) test data matrix.
  - `d`: Reduced dimension.
  - `B1`: Number of blocks of matrices.
  - `B2`: Number of blocks.
- **Outputs**:
  - `y_hat`: Estimated classes.

### 6. `RPE.R`
- **Key Function**: `RPE`
- **Inputs**:
  - `type`: Type of random matrix distribution (e.g., `'StdNormal'` for standard normal).
  - `X_train`: \( n \times p \) training data matrix.
  - `y_train`: Training classes.
  - `X_test`: \( m \times p \) test data matrix.
  - `d`: Reduced dimension.
- **Outputs**:
  - `y_hat`: Estimated classes.

## Scripts for Model Evaluation
To evaluate the misclassification probabilities for each method and model, the following scripts are used:

### 1. `Model.R`
- **Inputs**:
  - `type`: Scheme/Model type (e.g., `'model1'`, `'model2'`, `'model3'`, `'model4'`).
  - `p`: Data dimension.
- **Outputs**:
  - `mu1`, `mu2`: Mean parameters for the Gaussian distribution for class 1 and class 2.
  - `sigma1`, `sigma2`: Covariance matrices for the Gaussian distribution for class 1 and class 2.
  - `Omega1`, `Omega2`: Inverse of the covariance matrices for class 1 and class 2.
  - `logdet`: `logdet(sigma2) - logdet(sigma1)`.

### 2. `generateData_Bayes.R`
- **Inputs**:
  - Same as `Model.R`, plus:
  - `n1`, `n2`: Training observations for class 1 and class 2.
  - `m1`, `m2`: Test observations for class 1 and class 2.
- **Outputs**:
  - `X1`: \( (n1 + m1) \times p \) data matrix from class 1.
  - `X2`: \( (n2 + m2) \times p \) data matrix from class 2.
  - `res`: Bayes decision rule.

### 3. `Experiment.R`
- **Inputs**:
  - `model`: Schemes mentioned in the paper (e.g., `'model1'`, `'model2'`, `'model3'`, `'model4'`).
  - `method`: Methods mentioned in the paper (e.g., `'HDDA'`).
  - `p_all`: Vector of dimensions.
  - `iter`: Number of times data generation is executed.
  - `n1`, `n2`: Training observations for class 1 and class 2.
  - `m1`, `m2`: Test observations for class 1 and class 2.
  - `file_name`: File where outputs are stored (e.g., `'result.RData'`).
- **Outputs**:
  - `MP`: List of misclassification proportions of dimension `(iter, (n2 + m2), p_all)`.
  - `Ti`: List of execution times of dimension `(iter, (n2 + m2), p_all)`.

### 4. `Main_Parallel.R`
This script runs the full evaluation in a parallel computing environment.
