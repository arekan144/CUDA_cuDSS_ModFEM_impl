#pragma once
#ifndef SPARSESTRUCTURES
#define SPARSESTRUCTURES
#include <fstream>
#include <thread>
#include <future>
namespace SparseStructures
{
	/*
	* Simple wrapper for CSR - compresed sparse row format
	* sizeof(int) = 4
	*/
	class CSR {
	public:
		// Default constructor creates empty matrix, ready for asigment
		CSR() :a_csr(nullptr), row_ptr(nullptr), col_ind(nullptr), n(0u), nnz(0u) {};
		// Constructor for quick initalization when all values are already in memory
		CSR(double* a_csr, int* row_ptr, int* col_ind, unsigned int n, unsigned int nnz) :
			a_csr(a_csr), row_ptr(row_ptr), col_ind(col_ind), n(n), nnz(nnz) {};
		// Deep copy constructor
		CSR(const CSR& _CSR);
		// Deep assigment operator
		SparseStructures::CSR& operator=(const SparseStructures::CSR& _CSR);
		bool operator==(const SparseStructures::CSR& _CSR);
		// Reads file in special, 4 lines format N, NNZ, ROW, COL, VAL 
		static void readModFEMcrsMatrixFromFile(CSR& _CSR, std::ifstream& matrix_file);
		// Naive binary save. No compression. N, NNZ, ROW, COL, VAL
		static void saveModFEMcsrMatrixToBinary(CSR& _CSR, std::ofstream& matrix_file);
		// Reads saved non compressed binary files. In order: sizeof(int) N, sizeof(int) NNZ, sizeof(int)*(N+1) row pointers, sizeof(int)*NNZ column indices, sizeof(double)*NNZ non zero values
		static void readModFEMcsrMatrixFromBinary(CSR& _CSR, std::ifstream& matrix_file);

		double* getAcsr()
		{
			return a_csr;
		}
		int* getRowOff()
		{
			return row_ptr;
		}
		int* getColInd()
		{
			return col_ind;
		}
		int getN()
		{
			return n;
		}
		int getNNZ()
		{
			return nnz;
		}

		~CSR();
#ifdef _DEBUG
		static void print(const SparseStructures::CSR& csr, int max_nmb = 10);
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
#endif