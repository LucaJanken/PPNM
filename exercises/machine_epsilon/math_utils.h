#ifndef MATH_UTILS_H
#define MATH_UTILS_H

// Function to determine the maximum representable integer
int find_max_int();

// Function to determine the minimum representable integer
int find_min_int();

// Function to calculate double machine epsilon
double double_epsilon();

// Function to calculate float machine epsilon
float float_epsilon();

// Function to sum over epsilon/2 using method A
double sum_method_A(int n, double epsilon);

// Function to sum over epsilon/2 using method B
double sum_method_B(int n, double epsilon);

// Function to compare two double precision numbers with a given accuracy
bool approx(double a, double b, double acc = 1e-9, double eps = 1e-9);

#endif