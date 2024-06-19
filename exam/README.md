# Eigenvalues with Rayleigh Quotient and Locally Optimized Gradient Descent

## Introduction
This project implements a method to find the first and second energy eigenvalues of a symmetric matrix by minimizing the Rayleigh quotient using Newton's method with an analytic gradient and optimized line search. The hydrogen atom's ground state energy and the second smallest eigenvalue are computed using the grid representation Hamiltonian (as was also used in the EVD exercise).

## Code Overview

### Finding Eigenvalues
The code consists of several components organized into different files:

1. **Rayleigh Quotient Calculation**: Computes the Rayleigh quotient and its gradient for a given symmetric matrix and vector.
2. **Newton's Method**: Utilizes an analytic gradient and optimized line search to iteratively minimize the Rayleigh quotient and find the smallest eigenvalue and corresponding eigenvector.
3. **Deflation Procedure**: After finding the smallest eigenvalue, modifies the matrix to exclude this eigenvalue, allowing the algorithm to find the next smallest eigenvalue.
4. **Hydrogen Hamiltonian**: Uses a grid representation Hamiltonian previously used in the eigenvalue decomposition (EVD) homework exercise.
5. **Scaling Investigation**: Analyzes the time complexity of the algorithm by measuring the time taken to find eigenvalues for matrices of increasing size. Includes non-convergence failsafes and retries to handle ill-behaved matrices.

### How the Code Works

#### Finding the First Eigenvalue
1. **Rayleigh Quotient Calculation**: Computes the Rayleigh quotient \( R(H, v) = \frac{v^T H v}{v^T v} \) and its gradient \( \nabla R(H, v) \).
2. **Newton's Method**: Iteratively minimizes the Rayleigh quotient using the gradient and an optimized line search to determine the step size.
3. **Normalization**: Ensures the eigenvector remains normalized during the iterations.

#### Finding the Second Eigenvalue
1. **Deflation Procedure**: Projects the original matrix to remove the influence of the first eigenvector, creating a new matrix with the first eigenvalue removed.
2. **Rayleigh Quotient Calculation and Newton's Method**: Applies the same method to the new matrix to find the second smallest eigenvalue and corresponding eigenvector.

### Hydrogen Hamiltonian
The Hamiltonian for the hydrogen atom is generated using a grid representation. This approach was previously used in the EVD homework exercise. The smallest and second smallest eigenvalues are found using the methods described above.

### Scaling Investigation
The scaling investigation measures the time taken to compute the smallest eigenvalue for matrices of varying sizes. The results are written to `scaling_results.txt` and plotted in `scaling.svg` with an \( x^2 \) fit for comparison. The implementation includes non-convergence failsafes and a maximum retries parameter (currently set to 2) to discard ill-behaved randomly generated symmetric matrices.

### Code Execution
The code performs the following steps:
1. Computes and outputs the smallest eigenvalue and corresponding eigenvector for a test matrix.
2. Investigates the scaling of the method and writes results to `scaling_results.txt`.
3. Computes the ground state energy of the hydrogen atom and verifies the eigenvector normalization.
4. Finds the second smallest eigenvalue and corresponding eigenvector for the hydrogen atom.

The step-by-step procedure and results are written to `Out.txt`.

### Execution Time and Adjustments
The code is currently set to take approximately 20 seconds to build, primarily due to the scaling investigation and eigenenergy calculations. The execution time can be adjusted by:
- Reducing the upper bound for the scaling results.
- Increasing the grid spacing or reducing the maximum radius for the hydrogen Hamiltonian, which affects the accuracy of the eigenenergy approximation.

## Self-Evaluation
I would rate my efforts as 9/10. The exercise was completed as requested, but there is room for improvement. Specifically, a more robust algorithm that is less sensitive to high-dimensional matrices would have eliminated the need for failsafe procedures in the scaling investigation.
