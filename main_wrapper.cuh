#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include "matrixes/SparseStructures.hpp"
#include <cudss.h>
// for more debug info set CUDSS_LOG_LEVEL=5 in debug env
extern int doDecomposition(SparseStructures::CSR& matrix, double* b,
    double** x, cudssMatrixType_t mtype = CUDSS_MTYPE_GENERAL,
    cudssMatrixViewType_t mview = CUDSS_MVIEW_FULL,
    cudssIndexBase_t base = CUDSS_BASE_ZERO);