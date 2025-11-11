#include <iostream>
#include <fstream>
#include "main_wrapper.cuh"
#include "matrixes/SparseStructures.hpp"
#include <random>



int main()
{
    std::ifstream matrix_file("./matrixes/modfem_crs_1210.txt");
    if (!matrix_file.is_open()) {
        std::cout << "error\n";
        return -1;
    }

    SparseStructures::CSR test_csr = SparseStructures::CSR::readModFEMcrsMatrix(matrix_file);
    // SparseStructures::CSR::print(test_csr);

    matrix_file.close();

    double* b = new double[test_csr.getN()];

    for (unsigned int i = 0; i < test_csr.getN(); i++)
    {
        b[i] = 0.23;
    }

    double* x;

    doDecomposition(test_csr, b, &x);

    for (unsigned i = 0u; i < test_csr.getN(); i++)
        std::cout << x[i] << " ";
    std::cout << std::endl;
    delete[] b;
    delete[] x;

    return 0;
}
