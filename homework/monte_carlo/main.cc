#include <fstream>
#include <iomanip>
#include "montecarlo.h"

// Function to check if a point is inside a unit circle
double circle_area(vector& point) {
    // point[0] = x, point[1] = y
    if(point[0]*point[0] + point[1]*point[1] <= 1.0) return 1.0;
    else return 0.0;
}

// Function to check if a points is inside an ellipse
double ellipse_area(vector& point, double a) {
    // point[0] = x, point[1] = y
    if(point[0]*point[0] / (a*a) + point[1]*point[1] <= 1.0) return 1.0;
    else return 0.0;
}

double singularIntegrand(vector& point) {
    double x = point[0], y = point[1], z = point[2];
    double cosProd = std::cos(x) * std::cos(y) * std::cos(z);
    if (cosProd == 1) return 0; // Avoid division by zero
    return 1.0 / (1.0 - cosProd);
}

int main() {
    // Open .csv files for writing
    std::ofstream circleFile("circle_data.csv");
    std::ofstream ellipseFile("ellipse_data.csv");
    std::ofstream circleFileQuasi("circle_data_quasi.csv");
    std::ofstream ellipseFileQuasi("ellipse_data_quasi.csv");

    int nmin = 1500; // Threshold for stratified sampling

    // Define constants and bounds
    const double pi = 3.141592653589793;
    const double a = 2.0; // semi-major axis for ellipse
    const double trueAreaCircle = pi;
    const double trueAreaEllipse = pi * a;

    vector llCircle(2), urCircle(2); 
    llCircle[0] = -1.0; llCircle[1] = -1.0;
    urCircle[0] = 1.0; urCircle[1] = 1.0;

    vector llEllipse(2), urEllipse(2);
    llEllipse[0] = -a; llEllipse[1] = -1;
    urEllipse[0] = a; urEllipse[1] = 1;

    for (int N = 100; N <= 25000; N += 100) {
        // Calculate for circle
        auto resultCircle = plainmc(circle_area, llCircle, urCircle, N);
        double actualErrorCircle = std::abs(resultCircle.first - trueAreaCircle);

        // Calculate for ellipse
        auto resultEllipse = plainmc(std::bind(ellipse_area, std::placeholders::_1, a), llEllipse, urEllipse, N);
        double actualErrorEllipse = std::abs(resultEllipse.first - trueAreaEllipse);

        // Write to files
        circleFile << N << " " << resultCircle.first << " " << resultCircle.second << " " << actualErrorCircle << "\n";
        ellipseFile << N << " " << resultEllipse.first << " " << resultEllipse.second << " " << actualErrorEllipse << "\n";

        // Compute with quasi_mc for circle
        auto resultCircleQuasi = quasi_mc(circle_area, llCircle, urCircle, N);
        double actualErrorCircleQuasi = std::abs(resultCircleQuasi.first - trueAreaCircle);

        // Compute with quasi_mc for ellipse
        auto resultEllipseQuasi = quasi_mc(std::bind(ellipse_area, std::placeholders::_1, a), llEllipse, urEllipse, N);
        double actualErrorEllipseQuasi = std::abs(resultEllipseQuasi.first - trueAreaEllipse);

        // Write new results to quasi-random files
        circleFileQuasi << N << " " << resultCircleQuasi.first << " " << resultCircleQuasi.second << " " << actualErrorCircleQuasi << "\n";
        ellipseFileQuasi << N << " " << resultEllipseQuasi.first << " " << resultEllipseQuasi.second << " " << actualErrorEllipseQuasi << "\n";


        if (N == 25000) {
            std::cout << "Estimated area of circle: " << resultCircle.first << ", Estimated error: " << resultCircle.second << ", Actual error: " << actualErrorCircle << std::endl;
            std::cout << "Estimated area of ellipse: " << resultEllipse.first << ", Estimated error: " << resultEllipse.second << ", Actual error: " << actualErrorEllipse << std::endl;
            std::cout << "Estimated area of circle (quasi): " << resultCircleQuasi.first << ", Estimated error: " << resultCircleQuasi.second << ", Actual error: " << actualErrorCircleQuasi << std::endl;
            std::cout << "Estimated area of ellipse (quasi): " << resultEllipseQuasi.first << ", Estimated error: " << resultEllipseQuasi.second << ", Actual error: " << actualErrorEllipseQuasi << std::endl;
        }
    }


    // Close files
    circleFile.close();
    ellipseFile.close();
    circleFileQuasi.close();
    ellipseFileQuasi.close();

    int N = 25000; 

    // Calculate circle and ellipse areas with stratified sampling
    auto resultCircleStrat = stratified_mc(circle_area, llCircle, urCircle, N, nmin);
    double actualErrorCircleStrat = std::abs(resultCircleStrat.first - trueAreaCircle);
    auto resultEllipseStrat = stratified_mc(std::bind(ellipse_area, std::placeholders::_1, a), llEllipse, urEllipse, N, nmin);
    double actualErrorEllipseStrat = std::abs(resultEllipseStrat.first - trueAreaEllipse);
    // Calculate the difficult singular integral
    vector lowerBounds(3), upperBounds(3);
    lowerBounds[0] = lowerBounds[1] = lowerBounds[2] = 0.0;
    upperBounds[0] = upperBounds[1] = upperBounds[2] = pi;


    auto singularResult = plainmc(singularIntegrand, lowerBounds, upperBounds, N);
    auto singularResultQuasi = quasi_mc(singularIntegrand, lowerBounds, upperBounds, N);
    auto singularResultStrat = stratified_mc(singularIntegrand, lowerBounds, upperBounds, N, nmin);
    double scale_factor = 1.0 / (pi * pi * pi); // Scaling factor for the integral
    double actualErrorSingular = std::abs(singularResult.first * scale_factor - 1.39320392968567);
    double actualErrorSingularQuasi = std::abs(singularResultQuasi.first * scale_factor - 1.39320392968567);
    double actualErrorSingularStrat = std::abs(singularResultStrat.first * scale_factor - 1.39320392968567);

    std::cout << "Estimated area of circle (strat): " << resultCircleStrat.first << ", Estimated error: " << resultCircleStrat.second << ", Actual error: " << actualErrorCircleStrat << std::endl;
    std::cout << "Estimated area of ellipse (strat): " << resultEllipseStrat.first << ", Estimated error: " << resultEllipseStrat.second << ", Actual error: " << actualErrorEllipseStrat << std::endl;
    std::cout << "Estimated integral (singular): " << singularResult.first * scale_factor << ", Estimated error: " << singularResult.second * scale_factor << ", Actual error: " << actualErrorSingular << std::endl;
    std::cout << "Estimated integral (singular, quasi): " << singularResultQuasi.first * scale_factor << ", Estimated error: " << singularResultQuasi.second * scale_factor << ", Actual error: " << actualErrorSingularQuasi << std::endl;
    std::cout << "Estimated integral (singular, strat): " << singularResultStrat.first * scale_factor << ", Estimated error: " << singularResultStrat.second * scale_factor << ", Actual error: " << actualErrorSingularStrat << std::endl;

    return 0;
}
