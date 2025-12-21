#include "SparseStructures.h"
#include <string>

void SparseStructures::CSR::copyMatrix(const SparseStructures::CSR* inpt, SparseStructures::CSR* dest)
{
    if (inpt == dest) return;

    dest->n = inpt->n;
    dest->nnz = inpt->nnz;
    if (inpt->n == 0 || inpt->row_ptr == nullptr)
    {
        dest->col_ind = dest->row_ptr = nullptr;
        dest->a_csr = nullptr;
    }
    else {
        if (dest->a_csr != nullptr)
        {
            delete[] dest->a_csr;
            delete[] dest->row_ptr;
            delete[] dest->col_ind;
        }

        dest->a_csr = new double[static_cast<size_t>(inpt->nnz)];
        dest->row_ptr = new int[static_cast<size_t>(inpt->n + 1)];
        dest->col_ind = new int[static_cast<size_t>(inpt->nnz)];

        for (auto i = 0; i < inpt->n + 1; i++) {
            dest->row_ptr[i] = inpt->row_ptr[i];
        }

        for (auto i = 0; i < inpt->nnz; i++) {
            dest->col_ind[i] = inpt->col_ind[i];
            dest->a_csr[i] = inpt->a_csr[i];
        }
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
#ifdef _DEBUG
#include <iostream>
#endif
bool SparseStructures::CSR::operator==(const SparseStructures::CSR& _CSR)
{
    if (_CSR.nnz != this->nnz)
        return false;
    for (auto i = 0; i < _CSR.nnz; i++)
        if (std::abs(_CSR.a_csr[i] - this->a_csr[i]) > 1e-15) {
#ifdef _DEBUG
            std::cout << "Data not matching on: "<< i << std::endl;
#endif
            return false;
        }
    return true;
}

static void assignValuesInt(std::string& text_data_from_matrix, int maxN,
    unsigned long long prev_find, int* data)
{
    auto find = 0llu;
    auto n_str = 0llu;
    for (auto i = 0; i < maxN; i++) {
        find = text_data_from_matrix.find(' ', prev_find + 1);
        n_str = find - prev_find;
        data[i] = std::stoi(text_data_from_matrix.substr(prev_find + 1, n_str));
        prev_find = find;
    }
    

}

static void assignValuesDouble(std::string& text_data_from_matrix, int startN, int maxN, int jump,
    unsigned long long prev_find, double* data)
{
    auto find = 0llu;
    auto n_str = 0llu;
    for (auto i = startN; i < maxN; i += jump) {
        if(i != 0)
            for (auto j = 1; j < jump; j++) 
                prev_find = find = text_data_from_matrix.find(' ', prev_find + 1);
        find = text_data_from_matrix.find(' ', prev_find + 1);
        n_str = find - prev_find;
        data[i] = std::stod(text_data_from_matrix.substr(prev_find + 1, n_str));

        prev_find = find;
    }

}
#define THREADED
void SparseStructures::CSR::readModFEMcrsMatrixFromFile(CSR& _CSR, std::ifstream& matrix_file)
{
  matrix_file.seekg(0, matrix_file.end);
  std::string text_data_from_matrix(matrix_file.tellg(), 0);
  matrix_file.seekg(0);
  matrix_file.read(const_cast<char*>(text_data_from_matrix.data()), text_data_from_matrix.size());

  if (_CSR.nnz != 0) {
      delete[] _CSR.col_ind;
      delete[] _CSR.row_ptr;
      delete[] _CSR.a_csr;
  }

  unsigned long long prev_find = text_data_from_matrix.find('\n');
  _CSR.n = std::stoi(text_data_from_matrix.substr(0, prev_find));
  unsigned long long find = text_data_from_matrix.find('\n', prev_find + 1);
  unsigned long long n_str = find - prev_find;
  _CSR.nnz = std::stoi(text_data_from_matrix.substr(prev_find+1, n_str));
  prev_find = find;
  if (_CSR.a_csr != nullptr)
  {
      delete[] _CSR.a_csr;
      delete[] _CSR.row_ptr;
      delete[] _CSR.col_ind;
  }
  _CSR.a_csr = new double[_CSR.nnz] {0.};
  _CSR.row_ptr = new int[(_CSR.n + 1)];
  _CSR.col_ind = new int[_CSR.nnz];

#ifdef THREADED
  unsigned long long n_col_start = text_data_from_matrix.find('\n', prev_find + 1);
  unsigned long long n_a_csr_start = text_data_from_matrix.find('\n', n_col_start + 1);

  std::thread t3(assignValuesDouble, text_data_from_matrix, 0, _CSR.nnz, 1, n_a_csr_start, _CSR.a_csr);
  std::thread t1(assignValuesInt, text_data_from_matrix, _CSR.n + 1, prev_find, _CSR.row_ptr);
  std::thread t2(assignValuesInt, text_data_from_matrix, _CSR.nnz, n_col_start, _CSR.col_ind);

  t3.join();
  t2.join();
  t1.join();

  
#else
  //prev_find = last \n after nnz
  for (int i = 0; i < _CSR.n + 1; i++) {
      find = text_data_from_matrix.find(' ', prev_find + 1);
      n_str = find - prev_find;
      _CSR.row_ptr[i] = std::stoi(text_data_from_matrix.substr(prev_find + 1, n_str));
      prev_find = find;
  }

  //prev_find = n_col_start;
  //prev_find = last \n after last row_ptr
  for (int i = 0; i < _CSR.nnz; i++) {
      find = text_data_from_matrix.find(' ', prev_find + 1);
      n_str = find - prev_find;
      _CSR.col_ind[i] = std::stoi(text_data_from_matrix.substr(prev_find + 1, n_str));
      prev_find = find;
  }

  //prev_find = n_a_csr_start;
  //prev_find = last \n after last col_ind
  for (int i = 0; i < _CSR.nnz; i++) {
      find = text_data_from_matrix.find(' ', prev_find + 1);
      n_str = find - prev_find;
      _CSR.a_csr[i] = std::stod(text_data_from_matrix.substr(prev_find + 1, n_str));
      prev_find = find;
  }
#endif
  /*
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
  */

  return;
}

void SparseStructures::CSR::saveModFEMcsrMatrixToBinary(CSR& _CSR, std::ofstream& matrix_file)
{
    
    matrix_file.write(reinterpret_cast<char*>(&_CSR.n), sizeof(int));
    matrix_file.write(reinterpret_cast<char*>(&_CSR.nnz), sizeof(int));
    matrix_file.write(reinterpret_cast<char*>(_CSR.row_ptr), static_cast<size_t>(_CSR.n + 1) * sizeof(int));
    matrix_file.write(reinterpret_cast<char*>(_CSR.col_ind), static_cast<size_t>(_CSR.nnz) * sizeof(int));
    matrix_file.write(reinterpret_cast<char*>(_CSR.a_csr), static_cast<size_t>(_CSR.nnz) * sizeof(double));

}

void SparseStructures::CSR::readModFEMcsrMatrixFromBinary(CSR& _CSR, std::ifstream& matrix_file)
{
    if (_CSR.nnz != 0) {
        delete[] _CSR.col_ind;
        delete[] _CSR.row_ptr;
        delete[] _CSR.a_csr;
    }

    matrix_file.read(reinterpret_cast<char*>(&_CSR.n), sizeof(int));
    matrix_file.read(reinterpret_cast<char*>(&_CSR.nnz), sizeof(int));

    _CSR.a_csr = new double[static_cast<size_t>(_CSR.nnz)];
    _CSR.row_ptr = new int[static_cast<size_t>(_CSR.n + 1)];
    _CSR.col_ind = new int[static_cast<size_t>(_CSR.nnz)];

    matrix_file.read(reinterpret_cast<char*>(_CSR.row_ptr), static_cast<size_t>(_CSR.n + 1) * sizeof(int));
    matrix_file.read(reinterpret_cast<char*>(_CSR.col_ind), static_cast<size_t>(_CSR.nnz) * sizeof(int));
    matrix_file.read(reinterpret_cast<char*>(_CSR.a_csr), static_cast<size_t>(_CSR.nnz) * sizeof(double));

}

SparseStructures::CSR::~CSR()
{
    delete[] col_ind;
    delete[] row_ptr;
    delete[] a_csr;
}

#ifdef _DEBUG
#include <iostream>
void SparseStructures::CSR::print(const SparseStructures::CSR& _CSR, int max_nmb)
{
    std::cout << "DEBUG print\n";
    if (_CSR.n == 0) {
        std::cout << "EMPTY MATRIX\n";
        return;
    }
    std::cout << _CSR.n << " " << _CSR.nnz;
    std::cout << std::endl;
    for (int i = 0u; i < _CSR.n + 1 && i < max_nmb; i++) {
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


