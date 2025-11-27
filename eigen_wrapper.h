#pragma once

#include "matrixes/SparseStructures.h"
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

extern void eigenDecompositon(SparseStructures::CSR& matrix, double* b, double** x);
