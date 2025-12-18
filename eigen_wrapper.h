#pragma once
#include "matrixes/SparseStructures.h"
extern int saveSolution(std::ofstream& file, double* x, const int N);
extern int loadSolution(std::ifstream& file, double** x, const int N);
extern int eigenDecompositon(SparseStructures::CSR& matrix, const double* b, double** x);
