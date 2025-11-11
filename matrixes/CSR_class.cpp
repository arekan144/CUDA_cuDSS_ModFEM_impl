#include "SparseStructures.hpp"
#include <string>

void SparseStructures::CSR::copyMatrix(const SparseStructures::CSR* inpt, SparseStructures::CSR* dest)
{
    if (inpt == dest) return;
//std::cout << inpt->n << " " << inpt->nnz;
    dest->n = inpt->n;
    dest->nnz = inpt->nnz;
    dest->a_csr = new double[inpt->nnz];
    dest->row_ptr = new int[(inpt->n + 1lu)];
    dest->col_ind = new int[inpt->nnz];
    for (unsigned int i = 0u; i < inpt->n + 1lu; i++) {
        dest->row_ptr[i] = inpt->row_ptr[i];
    }

    for (unsigned int i = 0u; i < inpt->nnz; i++) {
        dest->col_ind[i] = inpt->col_ind[i];
        dest->a_csr[i] = inpt->a_csr[i];
    }
}

SparseStructures::CSR::CSR(const CSR& _CSR)
{
    SparseStructures::CSR::copyMatrix(&_CSR, this);
}
SparseStructures::CSR& SparseStructures::CSR::operator=(const SparseStructures::CSR& _CSR)
{
    SparseStructures::CSR::copyMatrix(&_CSR, this);
    return *this;
}
SparseStructures::CSR SparseStructures::CSR::readModFEMcrsMatrix(std::ifstream& matrix_file)
{
  CSR _CSR;

  matrix_file >> _CSR.n;
  matrix_file >> _CSR.nnz;
  
  _CSR.a_csr = new double[_CSR.nnz ] {0.};
  _CSR.row_ptr = new int[(_CSR.n + 1u) ] {0};
  _CSR.col_ind = new int[_CSR.nnz] {0};
  
  for(unsigned int i = 0u; i < _CSR.n + 1u; i++){
      matrix_file >> _CSR.row_ptr[i];
  }
  
  for (unsigned int i = 0u; i < _CSR.nnz; i++) {
      matrix_file >> _CSR.col_ind[i];
  }

  for (unsigned int i = 0u; i < _CSR.nnz; i++) {
      matrix_file >> _CSR.a_csr[i];
  }
//SparseStructures::CSR::print(_CSR);
  return _CSR;
}

SparseStructures::CSR::~CSR()
{
    delete[] col_ind;
    delete[] row_ptr;
    delete[] a_csr;
}

#ifdef _DEBUG
#include <iostream>
void SparseStructures::CSR::print(const SparseStructures::CSR& _CSR, unsigned int max_nmb)
{
    std::cout << "DEBUG print\n";
    if (_CSR.n == 0) {
        std::cout << "EMPTY MATRIX\n";
        return;
    }
    std::cout << _CSR.n << " " << _CSR.nnz;
    std::cout << std::endl;
    for (int i = 0u; i < _CSR.n + 1lu && i < max_nmb; i++) {
        std::cout <<  _CSR.row_ptr[i] << " ";
    }
    std::cout << std::endl;
    for (int i = 0u; i < _CSR.nnz && i < max_nmb; i++) {
        std::cout << _CSR.col_ind[i] << " ";
    }
    std::cout << std::endl;
    for (int i = 0u; i < _CSR.nnz && i < max_nmb; i++) {
        std::cout << _CSR.a_csr[i] << " ";
    }
    std::cout << std::endl;
}
#endif


