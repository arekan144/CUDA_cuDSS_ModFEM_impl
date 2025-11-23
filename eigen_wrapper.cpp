#include "eigen_wrapper.h"
#include <Eigen/Sparse>
#include <Eigen/SparseLU.h>


using namespace Eigen;

int doDecompositonAndCompare(SparseStructures::CSR& matrix, double* b, double* x)
{
    SparseMatrix<double> A;
    
    SparseLU<SparseMatrix<double>, COLAMDOrdering<int> > solver;



    return 0;
}
