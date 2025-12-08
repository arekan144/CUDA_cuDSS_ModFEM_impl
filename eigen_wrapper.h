#pragma once

#include "matrixes/SparseStructures.h"
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <iostream>

extern int eigenDecompositon(SparseStructures::CSR& matrix, double* b, double** x);
