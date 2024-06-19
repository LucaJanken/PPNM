#include <iostream>
#include <iomanip>
#include <chrono>
#include "EigenSolver.h"
#include "HydrogenHamiltonian.h"
#include "ScalingInvestigation.h"

int main() {

    std::cout << "--------------------Test Eigenvalue Problem--------------------\n" << std::endl;

    // Define a symmetric matrix H and print it
    matrix H(3, 3);
    H.set(0, 0, 4); H.set(0, 1, 1); H.set(0, 2, 1);
    H.set(1, 0, 1); H.set(1, 1, 3); H.set(1, 2, 1);
    H.set(2, 0, 1); H.set(2, 1, 1); H.set(2, 2, 2);

    H.print("Test Matrix H:");

    std::cout << "\nStep 1: Finding the smallest eigenvalue and corresponding eigenvector of the matrix H using the Rayleigh Quotient method...\n";
    auto eigenpair = FindEigenPair(H, 1e-6, 5.0);

    std::cout << std::fixed << std::setprecision(9);

    std::cout << "\nStep 2: Output the smallest eigenvalue.\n";
    std::cout << "Smallest Eigenvalue: " << eigenpair.first << std::endl;

    std::cout << "\nStep 3: Output the corresponding eigenvector.\n";
    std::cout << "Corresponding Eigenvector: ";
    for (double vi : eigenpair.second) {
        std::cout << vi << " ";
    }
    std::cout << std::endl;

    std::cout << "\nStep 4: Verify and output the normalization of the eigenvector.\n";
    double norm = 0.0;
    for (double vi : eigenpair.second) {
        norm += vi * vi;
    }
    norm = std::sqrt(norm);
    std::cout << "Norm of the Eigenvector: " << norm << " (should be close to 1.0)" << std::endl;

    std::cout << "\nStep 5: Compute the product H * v, where v is the eigenvector.\n";
    vector Hv(3, 0.0);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            Hv[i] += H(i, j) * eigenpair.second[j];
        }
    }

    std::cout << "H * v: ";
    for (double hvi : Hv) {
        std::cout << hvi << " ";
    }
    std::cout << std::endl;

    std::cout << "\nStep 6: Compute the product lambda * v, where lambda is the smallest eigenvalue and v is the eigenvector.\n";
    std::cout << "lambda * v: ";
    for (double vi : eigenpair.second) {
        std::cout << eigenpair.first * vi << " ";
    }
    std::cout << std::endl;

    std::cout << "\nStep 7: Verification.\n";
    std::cout << "The values of 'H * v' and 'lambda * v' should be close to each other.\n";
    std::cout << "'H * v' represents the matrix-vector multiplication of H with the eigenvector v.\n";
    std::cout << "'lambda * v' represents the scalar multiplication of the eigenvalue with the eigenvector v.\n";
    std::cout << "If the values match closely, it confirms that the computed eigenpair is correct.\n";

    // Investigate scaling
    std::cout << "\n--------------------Investigate Scaling of the Method (See scaling_results.txt)--------------------\n" << std::endl;
    try {
        investigate_scaling();
    } catch (const std::exception &e) {
        std::cerr << "Scaling investigation stopped: " << e.what() << std::endl;
    }

    std::cout << "\n--------------------Calculate Ground State Energy of Hydrogen Atom--------------------\n" << std::endl;

    // Set parameters for the Hydrogen atom problem
    double dr = 0.3; // Grid spacing
    int rmax = 13;   // Maximum radius
    float hartree_to_ev = 27.2114; // Conversion factor from Hartree to eV

    // Generate the Hamiltonian matrix for Hydrogen
    std::cout << "\nStep 8: Generating the Hamiltonian matrix for the Hydrogen atom (rmax = " << rmax << ",  dr = " << dr << ")...\n";
    matrix Hyd = Hamiltonian(rmax, dr);
    
    // Find the smallest eigenvalue and corresponding eigenvector
    std::cout << "\nStep 9: Finding the smallest eigenvalue and corresponding eigenvector of the Hydrogen atom Hamiltonian...\n";
    auto eigenpairHyd = FindEigenPair(Hyd, 1e-6, 10.0);

    // Convert the smallest eigenvalue to eV
    double smallestEigenvalueHartree = eigenpairHyd.first;
    double smallestEigenvalueeV = smallestEigenvalueHartree * hartree_to_ev;

    // Output the results
    std::cout << "\nStep 10: Output the smallest eigenvalue of the Hydrogen atom in Hartree and eV.\n";
    std::cout << "Smallest Eigenvalue of Hydrogen: " << smallestEigenvalueHartree << " Hartree";
    std::cout << " (" << smallestEigenvalueeV << " eV, should be close to -13.6 eV)" << std::endl;

    std::cout << "\nStep 11: Verify the normalization of the eigenvector of the Hydrogen atom.\n";
    double normHyd = 0.0;
    for (double vi : eigenpairHyd.second) {
        normHyd += vi * vi;
    }
    normHyd = std::sqrt(normHyd);
    std::cout << "Norm of the Eigenvector of Hydrogen: " << normHyd << " (should be close to 1.0)" << std::endl;

    std::cout << "\n--------------------Calculate Second Smallest Eigenergy of Hydrogen Atom--------------------\n" << std::endl;

    // Find the second smallest eigenvalue and corresponding eigenvector for Hydrogen
    std::cout << "\nStep 12: Finding the second smallest eigenvalue and corresponding eigenvector for the Hydrogen atom Hamiltonian...\n";
    auto second_eigenpairHyd = FindSecondEigenPair(Hyd, eigenpairHyd.second, 1e-6, 10.0);

    // Convert the second smallest eigenvalue to eV
    double secondEigenvalueHartree = second_eigenpairHyd.first;
    double secondEigenvalueeV = secondEigenvalueHartree * hartree_to_ev;

    // Output the results for the second smallest eigenvalue
    std::cout << "Second Smallest Eigenvalue of Hydrogen: " << secondEigenvalueHartree << " Hartree";
    std::cout << " (" << secondEigenvalueeV << " eV, should be close to -3.4 eV)" << std::endl;

    // Verify the normalization of the eigenvector
    double normSecondHyd = 0.0;
    for (double vi : second_eigenpairHyd.second) {
        normSecondHyd += vi * vi;
    }
    normSecondHyd = std::sqrt(normSecondHyd);
    std::cout << "Norm of the Eigenvector of Hydrogen (Second Smallest Eigenvalue): " << normSecondHyd << " (should be close to 1.0)" << std::endl;

    return 0;
}
