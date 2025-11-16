#include <iostream>
#include <fstream>
#include <random>

#include "main_wrapper.cuh"
#include "matrixes/SparseStructures.hpp"
#include "ChronoTimer.h"

int main()
{
    std::ifstream matrix_file("./matrixes/modfem_crs_3745223.txt");
    if (!matrix_file.is_open()) {
        std::cout << "error\n";
        return -1;
    }

    ChronoTimer timer;

    SparseStructures::CSR test_csr;
timer.setStartTime();
    SparseStructures::CSR::readModFEMcrsMatrixFromFile(test_csr, matrix_file);
timer.setEndTime();

std::cout << "Time difference: "<< timer.getDiffInMS() << " [ms]" << std::endl;
std::cout << "Time difference: " << timer.getDiffInS() << " [s]" << std::endl;
   
    matrix_file.close();

#ifdef _DEBUG
    SparseStructures::CSR::print(test_csr,100);
#endif
//test
    return 0;

    double* b = new double[test_csr.getN()];

    for (int i = 0; i < test_csr.getN(); i++)
        b[i] = 1;

    double* x;

    doDecomposition(test_csr, b, &x);

    for (int i = 0; i < test_csr.getN(); i++)
        std::cout << x[i] << " ";
    std::cout << std::endl;
    delete[] b;
    delete[] x;

    return 0;
}
