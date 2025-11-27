#include <iostream>
#include <fstream>
#include <random>

#include "cudss_wrapper.cuh"
#include "eigen_wrapper.h"
#include "matrixes/SparseStructures.h"
#include "ChronoTimer.h"

int main() {

    std::ifstream matrix_file("./matrixes/modfem_crs_1210.txt");
    if (!matrix_file.is_open()) {
        std::cout << "error\n";
        return -1;
    }
    
    SparseStructures::CSR test_csr;
    SparseStructures::CSR::readModFEMcrsMatrixFromFile(test_csr, matrix_file);

    matrix_file.close();

#ifdef _DEBUG
    SparseStructures::CSR::print(test_csr,100);
#endif

    double* b = new double[test_csr.getN()];

    for (auto i = 0u; i < test_csr.getN(); i++)
        b[i] = 1.;

    double* x_1, *x_2;

    eigenDecompositon(test_csr, b, &x_2);

    size_t how_much = 10Ui64;
    ChronoTimer timer(6UI64 * how_much);
    for (auto i = 0Ui64; i < how_much ; i++)
    {
        cuDSSDecompositionWithSynchronization(timer, test_csr, b, &x_1, CUDSS_MTYPE_GENERAL, CUDSS_MVIEW_FULL, CUDSS_BASE_ZERO);
        auto tolerance = 1.e-15;
        int times = 0;
        for (auto i = 0u; i < test_csr.getN(); i++)
            if (std::abs(x_1[i] - x_2[i]) < tolerance)
                times++;

        std::cout << "For attempt " << i << "\n";
        std::cout << "Inside given tolerance: " << tolerance << " " << times << "\npercent " << static_cast<double>(times) / static_cast<double>(test_csr.getN()) << " \n";
    }
    
    std::ofstream timefile("testtime.txt");
    if (timefile.is_open()) {
        timer.saveTimesToFile(timefile,6);
        timefile.close();
    }
    delete[] b;
    delete[] x_1;
    delete[] x_2;

    return 0;
}
