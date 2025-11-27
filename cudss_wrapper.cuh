#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include "matrixes/SparseStructures.h"
#include <cudss.h>
#include "ChronoTimer.h"
// for more debug info set CUDSS_LOG_LEVEL=5 in debug env

/*cuDSSDecompositionWithSynchronization makes synchronization calls (cudaDeviceSynchronize)
  individual testing of:
  * CUDA[init,aloc,copy] + cuDSS[init], 
  * cuDSS[creation of wrappers of sparse and dense matrixes]
  * cuDSS[execute: analisys]
  * cuDSS[execute: factorization]
  * cuDSS[execute: solve]
  * CUDA[copy]
  Recomended timer size = above count (6)
*/ 
extern int cuDSSDecompositionWithSynchronization(ChronoTimer& timer, SparseStructures::CSR& matrix,
    double* b, double** x,
    cudssMatrixType_t mtype = CUDSS_MTYPE_GENERAL,
    cudssMatrixViewType_t mview = CUDSS_MVIEW_FULL, cudssIndexBase_t base = CUDSS_BASE_ZERO);