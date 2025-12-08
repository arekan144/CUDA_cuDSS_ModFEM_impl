#include "eigen_wrapper.h"

using namespace Eigen;
int eigenDecompositon(SparseStructures::CSR& matrix, double* b, double** x)
{
    try {
        // Copy CSR to Eigen SparseMatrix, utilizing built-in Map class
        Eigen::Map<Eigen::SparseMatrix<double, ColMajor>> mat_map(
            matrix.getN(), matrix.getN(), matrix.getNNZ(),
            matrix.getRowOff(), matrix.getColInd(), matrix.getAcsr());
        // Eval tranfroms given matrix to SparseMatrix format
        SparseMatrix<double, ColMajor> A = mat_map.eval();
        // "SparseMatrixes should be in a compressed and column-major form" - excerpt from oficial documentation
        if (!A.isCompressed())
            A.makeCompressed();
        // One row (column major) matrixes
        VectorXd x_e(matrix.getN()), b_e(matrix.getN());

        for (auto i = 0; i < matrix.getN(); i++)
            b_e[i] = b[i];

        // solver initialization and computation with given parameters
        SparseLU<SparseMatrix<double>, COLAMDOrdering<int> > solver;
        solver.analyzePattern(A);
        solver.factorize(A);
        x_e = solver.solve(b_e);

        *x = new double[matrix.getN()];

        for (auto i = 0; i < matrix.getN(); i++)
            (*x)[i] = x_e[i];
    }
    catch (std::exception exc) {
        std::cout << exc.what() << std::endl;
        return 1;
    }
    return 0;
}
