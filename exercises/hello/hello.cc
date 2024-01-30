#include <iostream>                                             // Library for in- and output operations
#include <string>                                               // Library for string manipulation

int main() {                                                    // Convention in C++ for every program to start with 'main'
    std::string username = "LucaJanken";                        // 'std' Standard namespace.
    std::cout << "Hello, " << username << "!" << std::endl;     // 'std::cout' Output stream object (for writing to console). '<<' to concatenate. 'std::endl' newline character
    return 0;                                                   // Used to indicate a successful program execution. typically '0' for success and non-zero values for errors.
}