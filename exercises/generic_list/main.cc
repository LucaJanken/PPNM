#include "genlist.h"
#include <iostream>

int main() {
    genlist<double> list;
    double x;

    while (std::cin >> x) {
        std::cerr << "x=" << x << std::endl;
        list.add(x);
    }

    for (size_t i = 0; i < list.size; i++) {
        std::cout << "list[" << i << "]=" << list[i] << std::endl;
    }

    genlist<genlist<double>> dlist;
    dlist.add(list);
    genlist<double> d1;
    d1.add(5.4e-2);
    dlist.add(d1);
    d1.add(6.7e8);
    dlist.add(d1);
    d1.add(0.99);
    dlist.add(d1);

    std::cout << "dlist=" << std::endl;
    for (size_t i = 0; i < dlist.size; i++) {
        for (size_t j = 0; j < dlist[i].size; j++) {
            std::cout << dlist[i][j] << " ";
        }
        std::cout << std::endl;
    }
}
