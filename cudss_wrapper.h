#include "matrixes/SparseStructures.h"
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
    double* b, double** x, unsigned short matrix_type = 0, unsigned short view_type = 0, unsigned index_base = 0);

/*
typedef enum cudssMatrixType_t {
    CUDSS_MTYPE_GENERAL, 0
    CUDSS_MTYPE_SYMMETRIC, 1
    CUDSS_MTYPE_HERMITIAN, 2
    CUDSS_MTYPE_SPD, 3
    CUDSS_MTYPE_HPD 4
} cudssMatrixType_t;

typedef enum cudssMatrixViewType_t {
    CUDSS_MVIEW_FULL, 0
    CUDSS_MVIEW_LOWER, 1
    CUDSS_MVIEW_UPPER 2
} cudssMatrixViewType_t;

typedef enum cudssIndexBase_t {
    CUDSS_BASE_ZERO, 0
    CUDSS_BASE_ONE  1
*/
