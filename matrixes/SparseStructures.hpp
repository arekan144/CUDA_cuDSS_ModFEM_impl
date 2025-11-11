#pragma once
#include "fstream"
namespace SparseStructures
{
	class CSR {
	public:
		CSR() :a_csr(nullptr), row_ptr(nullptr), col_ind(nullptr), n(0u), nnz(0u) {};
		CSR(double* a_csr, int* row_ptr, int* col_ind, unsigned int n, unsigned int nnz) :
			a_csr(a_csr), row_ptr(row_ptr), col_ind(col_ind), n(n), nnz(nnz) {
		};
		CSR(const CSR& _CSR);
		SparseStructures::CSR& operator=(const SparseStructures::CSR& _CSR);
		static SparseStructures::CSR readModFEMcrsMatrix(std::ifstream& matrix_file);
		double* SparseStructures::CSR::getAcsr()
		{
			return a_csr;
		}
		int* SparseStructures::CSR::getRowPtr()
		{
			return row_ptr;
		}
		int* SparseStructures::CSR::getColInd()
		{
			return col_ind;
		}
		unsigned int SparseStructures::CSR::getN()
		{
			return n;
		}
		unsigned int SparseStructures::CSR::getNNZ()
		{
			return nnz;
		}
		~CSR();
#ifdef _DEBUG
		static void print(const SparseStructures::CSR& csr, unsigned int max_nmb = 10);
#else
		static void print(const SparseStructures::CSR& csr) = delete;
#endif //_DEBUG
	private:
		static void copyMatrix(const SparseStructures::CSR* inpt, SparseStructures::CSR* dest);

		double* a_csr;		// main value table
		int* row_ptr;		// row pointer table
		int* col_ind;		// column indicies table
		int n;				// real(non-compressed) size of the matrix
		int nnz;			// number of non zero values
	};
}

