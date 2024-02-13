#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>

void processNumber(double x) {
    std::cout << "x = " << x << ", sin(x) = " << sin(x) << ", cos(x) = " << cos(x) << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc > 1) {
        std::string infile, outfile;
        bool numbersFound = false;

        for (int i = 1; i < argc; ++i) {
            std::string arg = argv[i];
            size_t pos = arg.find(":");
            if (pos != std::string::npos) {
                std::string key = arg.substr(0, pos);
                std::string value = arg.substr(pos + 1);

                if (key == "-numbers") {
                    std::stringstream ss(value);
                    std::string number;
                    while (std::getline(ss, number, ',')) {
                        double x = std::stod(number);
                        processNumber(x);
                    }
                    numbersFound = true;
                } else if (key == "-input") {
                    infile = value;
                } else if (key == "-output") {
                    outfile = value;
                }
            }
        }

        if (!infile.empty() && !outfile.empty()) {
            std::ifstream instream(infile);
            std::ofstream outstream(outfile, std::ofstream::out);
            if (!instream.is_open() || !outstream.is_open()) {
                std::cerr << "Error opening file" << std::endl;
                return 1;
            }

            double x;
            while (instream >> x) {
                outstream << "x = " << x << ", sin(x) = " << sin(x) << ", cos(x) = " << cos(x) << std::endl;
            }

            instream.close();
            outstream.close();
            return 0;
        } else if (!numbersFound) {
            // Fallback to reading from standard input if no specific command-line arguments are provided.
            std::string line;
            while (std::getline(std::cin, line)) {
                std::istringstream iss(line);
                double x;
                while (iss >> x) {
                    processNumber(x);
                }
            }
        }
    } else {
        // Default to reading from standard input if no arguments are given.
        std::string line;
        while (std::getline(std::cin, line)) {
            std::istringstream iss(line);
            double x;
            while (iss >> x) {
                processNumber(x);
            }
        }
    }

    return 0;
}
