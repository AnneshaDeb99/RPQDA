# RPQDA

The following description of the functions refers to the paper titled *"Ultrahigh-dimensional Quadratic Discriminant Analysis Using Random Projections"*.

## Function Descriptions

### 1. `RPE.R`
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
To evaluate the misclassification probabilities for each method and scheme, the following scripts are used:

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
  - `model`: Schemes mentioned in the paper (e.g., `'scheme1'`, `'scheme2'`, `'scheme3'`, `'scheme4'`).
  - `method`: Methods mentioned in the paper (e.g., `'HDDA'`, `'AoYa'`, `'DA-QDA'`, `'IIS-SQDA'`, `'RPE-CS'`, `'RPE-SN'`, `'RPE-TP'`).
  - `p_all`: Vector of dimensions.
  - `iter`: Number of times data generation is executed.
  - `n1`, `n2`: Training observations for class 1 and class 2.
  - `m1`, `m2`: Test observations for class 1 and class 2.
  - `file_name`: File where outputs are stored (e.g., `'result.RData'`).
- **Outputs**:
  - `MP`: array of misclassification proportions of dimension `(methods, iter, p_all)`.
  - `Ti`: List of execution times of dimension `(methods, iter, p_all)`.

### 4. `Main_Parallel.R`
This script runs the full evaluation in a parallel computing environment.

## How to Run the Codes
1. In `Main_Parallel.R`, specify the desired `scheme`, `method`, `p_all`, `iter`, and sample sizes (`n1`, `n2`, `m1`, `m2`).
2. Set the output directory where the results will be stored.

### Additional Details
- `Experiment.R` computes the misclassification proportion (`MP`) and execution time (`Ti`) for each method (e.g., Bayes, RPE-CS) across different dimensions and iterations.
- Within `Experiment.R`:
  - `Model.R` calculates the population parameters (`mu1`, `mu2`, `sigma1`, `sigma2`) based on the chosen scheme.
  - Using these parameters, `generateData_Bayes.R` generates data (`X1`, `X2`) and assigns classes based on the Bayes decision rule (`res`).

