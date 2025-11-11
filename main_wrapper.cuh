#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include "matrixes/SparseStructures.hpp"
#include <cudss.h>
// for more debug info set CUDSS_LOG_LEVEL=5 in debug env
extern int doDecomposition(SparseStructures::CSR& matrix, double* b, double** x);