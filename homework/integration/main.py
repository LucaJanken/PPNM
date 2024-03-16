import numpy as np
from scipy.integrate import quad

# Define the functions
def f_inv_sqrt(x):
    return 1 / np.sqrt(x)

def f_ln_sqrt(x):
    return np.log(x) / np.sqrt(x)

def integrand(x):
    return np.exp(-x)

#gaussian
def f_gaussian(x):
    return np.exp(-x**2)

# Integrate 1/sqrt(x) from 0 to 1
result_inv_sqrt, error_inv_sqrt = quad(f_inv_sqrt, 0, 1)
evaluations_inv_sqrt = quad(f_inv_sqrt, 0, 1, full_output=1)[2]['neval']

# Integrate ln(x)/sqrt(x) from 0 to 1
result_ln_sqrt, error_ln_sqrt = quad(f_ln_sqrt, 0, 1)
evaluations_ln_sqrt = quad(f_ln_sqrt, 0, 1, full_output=1)[2]['neval']

# Integrate exp(-x) from 0 to infinity
result_exp, error_exp = quad(integrand, 0, np.inf)
evaluations_exp = quad(integrand, 0, np.inf, full_output=1)[2]['neval']

# Integrate exp(-x^2) from -inf to inf
result_gaussian, error_gaussian = quad(f_gaussian, -np.inf, np.inf)
evaluations_gaussian = quad(f_gaussian, -np.inf, np.inf, full_output=1)[2]['neval']

print(f"(Python) Integral of 1/sqrt(x) over [0, 1]: Result = {result_inv_sqrt}, Evaluations = {evaluations_inv_sqrt}")
print(f"(Python) Integral of ln(x)/sqrt(x) over [0, 1]: Result = {result_ln_sqrt}, Evaluations = {evaluations_ln_sqrt}")
print(f"(Python) Integral of exp(-x) over [0, inf]: Result = {result_exp}, Evaluations = {evaluations_exp}, Error = {error_exp}")
print(f"(Python) Integral of exp(-x^2) over [-inf, inf]: Result = {result_gaussian}, Evaluations = {evaluations_gaussian}, Error = {error_gaussian}")
print("The C++ implementation seemingly requires more evaluations than numpy for these infinite limits, could perhaps be improved using the Clenshaw–Curtis variable transformation again.")