#include "eigen_wrapper.h"
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <iostream>
#include <fstream>
using namespace Eigen;

int saveSolution(std::ofstream& file, double* x, const int N)
{
	if(N < 0) return 1;
    	file.write(reinterpret_cast<char*>(const_cast<int*>(&N)), sizeof(int));
	file.write(reinterpret_cast<char*>(x), static_cast<size_t>(N) * sizeof(double));
	return 0;
}
int loadSolution(std::ifstream& file, double** x, const int N)
{
	if(N < 0) return 1;
	int conf = 0;
	file.read(reinterpret_cast<char*>(&conf), sizeof(int));
	if(conf != N) return 2;
	(*x) = new double[static_cast<size_t>(N)];
	file.read(reinterpret_cast<char*>((*x)), static_cast<size_t>(N) * sizeof(double));
	return 0;
}	

int eigenDecompositon(SparseStructures::CSR& matrix, const double* b, double** x)
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

        *x = new double[static_cast<size_t>(matrix.getN())];

        for (auto i = 0; i < matrix.getN(); i++)
            (*x)[i] = x_e[i];
    }
    catch (std::exception exc) {
        std::cout << exc.what() << std::endl;
        return 1;
    }
    return 0;
}
