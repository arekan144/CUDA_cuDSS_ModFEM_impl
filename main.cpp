#include <iostream>
#include <fstream>
#include <random>
#include <string>

#include "matrixes/SparseStructures.h"
#include "cudss_wrapper.h"
#include "eigen_wrapper.h"

#include "ChronoTimer.h"
#include "ArgumentParser.h"

#ifndef HOW_MUCH
#define HOW_MUCH 10llu 
#endif

static int __saveSolution(std::ofstream& file, double* x, const int N)
{
	if (N < 0) return 1;
	file.write(reinterpret_cast<char*>(const_cast<int*>(&N)), sizeof(int));
	file.write(reinterpret_cast<char*>(x), static_cast<size_t>(N) * sizeof(double));
	return 0;
}
static int __loadSolution(std::ifstream& file, double** x, const int N)
{
	if (N < 0) return 1;
	int conf = 0;
	file.read(reinterpret_cast<char*>(&conf), sizeof(int));
	if (conf != N) return 2;
	(*x) = new double[static_cast<size_t>(N)];
	file.read(reinterpret_cast<char*>((*x)), static_cast<size_t>(N) * sizeof(double));
	return 0;
}

static int __loadMatrix(std::string& file_name, short file_type, SparseStructures::CSR& _CSR) {
	std::ifstream matrix_file;
	
	switch (file_type) {
	case 1:
		matrix_file.open(file_name, std::ifstream::binary);
		SparseStructures::CSR::readModFEMcsrMatrixFromBinary(_CSR, matrix_file);
		break;
	case 0:
		matrix_file.open(file_name);
		SparseStructures::CSR::readModFEMcrsMatrixFromFile(_CSR, matrix_file);
		break;
	default:
		std::cout << "Not implemented\n";
		return -2;
	}
	std::cout << "File opended. Computing.\n";
	matrix_file.close();
	return 0;
}

static void __compare(SparseStructures::CSR& _CSR, double* x_1, double* x_2) {
	// we are looking for smallest difference between both solutions, 
	// where 100% of all values fall inside given treshold 
	const size_t data_points = 12llu;
	double tolerance[data_points] = { 1.e-16 };
	int times[data_points] = { 0 };
	for (auto i = 1llu; i < data_points; i++) {
		tolerance[i] = tolerance[i - 1] * 10.;
		times[i] = 0;
	}
	for (auto i = 0; i < _CSR.getN(); i++)
	{
		double abs_x = std::abs(x_1[i] - x_2[i]);
		for (auto k = 0llu; k < data_points; k++)
			if (abs_x < tolerance[k])
				times[k]++;
	}
	
	double next_percent = static_cast<double>(times[0]) * 100. / static_cast<double>(_CSR.getN());
	size_t k = 0llu;
	do {
		std::cout << "Inside given tolerance: " << tolerance[k] << " " << times[k] << ", percent: " << next_percent << " % \n";
		if (next_percent >= 100.) break;
		next_percent = static_cast<double>(times[++k]) * 100. / static_cast<double>(_CSR.getN());
	} while (k != data_points);
}

static int __getSparcityPatterns(std::string& file_name, short file_type, short matrix_type, int max_res = 0) {
	SparseStructures::CSR test_csr;
	__loadMatrix(file_name, file_type, test_csr);
	double* b = new double[test_csr.getN()];

	for (auto i = 0; i < test_csr.getN(); i++)
		b[i] = 1.;
	double* x_1 = nullptr, * x_2 = nullptr;
	std::cout << max_res << "\n";
	
	if(max_res)
		for(short algorithm = 0; algorithm < 4; algorithm++)
			cuDSSOnlyAnalisysAndSpPattern(test_csr, b, &x_1, matrix_type, algorithm, max_res);
	else
		for (short algorithm = 0; algorithm < 4; algorithm++)
			cuDSSOnlyAnalisysAndSpPattern(test_csr, b, &x_1, matrix_type, algorithm);
	delete[] b;
	delete[] x_1;
	return 0;
}

static int __cudssAndEigen(std::string& file_name, short file_type, short nskip_eigen, short matrix_type, std::string solve_name)
{
	SparseStructures::CSR test_csr;
	__loadMatrix(file_name, file_type, test_csr);

#ifdef _DEBUG
	SparseStructures::CSR::print(test_csr, 100);
#endif

	double* b = new double[test_csr.getN()];

	for (auto i = 0; i < test_csr.getN(); i++)
		b[i] = 1.;

	double* x_1 = nullptr, * x_2 = nullptr;
	switch (nskip_eigen) {
		// do not skip
	case static_cast<short>(1):
	case static_cast<short>(2):
		std::cout << "Started Eigen. ";
		eigenDecompositon(test_csr, b, &x_2);
		std::cout << "Done.\n";
		break;
		// provided file with previously calculated solution
	case static_cast<short>(3):
	{
		std::cout << "Opening previous solution: " << solve_name << std::endl;

		std::ifstream solution_file(solve_name, std::ifstream::binary);

		if (!solution_file.is_open()) {
			std::cout << "Cannot open file " << solve_name << ".\n Skipping eigenDecomposition.\n";
			nskip_eigen = static_cast<short>(0);
			break;
		}
		int ret = __loadSolution(solution_file, &x_2, test_csr.getN());
		if (ret) {
			std::cout << "File reading failed.\n Skipping eigenDecomposition.\n";
			nskip_eigen = static_cast<short>(0);
		}
		solution_file.close();
	}
	break;
	// skipping decomposition using Eigen library
	case static_cast<short>(0):
	default:
		std::cout << "Skipped eigenDecomposition" << "\n";
		break;
	}

	// How many times decomposition should be called. HOW_MUCH is defined at the begining of this file, if not provided by command line.
	constexpr size_t how_much = HOW_MUCH;
	ChronoTimer timer(6llu * how_much);
	for (auto i = 0llu; i < how_much; i++) {
		cuDSSDecompositionWithSynchronization(timer, test_csr, b, &x_1, matrix_type, 0, 0);
		if (!nskip_eigen) {
			continue;
		}
		std::cout << "For attempt " << i << "\n";
		__compare(test_csr, x_1, x_2);
	}

	size_t last_slash = file_name.rfind('/', file_name.size() - 1) + 1;
	std::string time_name = file_name.substr(last_slash, file_name.rfind(".txt") - last_slash) + "_" + std::to_string(matrix_type) + "_test_time.txt";

	std::ofstream timefile(time_name);
	if (timefile.is_open()) {
		timer.saveTimesToFile(timefile, 6);
		timefile.close();
		std::cout << "Time stamps saved under file: " << time_name << "\n";
	}

	if (nskip_eigen == static_cast<short>(2)) {
		std::string value_name = file_name.substr(last_slash, file_name.rfind(".txt") - last_slash) + "_solved";
		std::cout << "Saving result to txt file under " << value_name << ". ";
		std::ofstream solution_file(value_name, std::ofstream::binary);
		if (solution_file.is_open()) {
			__saveSolution(solution_file, x_2, test_csr.getN());
			solution_file.close();
			std::cout << "Done.\n";
		}
		else {
			std::cout << "Cannot open file: " << value_name << '\n';
		}

	}

	delete[] b;
	delete[] x_1;
	delete[] x_2;
	return 0;
}

static int __cudssAndCusolve(std::string& file_name, short file_type, short nskip_cusolve, short matrix_type, std::string solve_name)
{

	return 0;
}

static int __saveToBinary(std::string file_name) {
	std::cout << "Saving as binary from file " << file_name << std::endl;
	std::ifstream load_file(file_name, std::ifstream::binary);
	if (!load_file.is_open()) {
		std::cout << "Can't open file\n";
		return -1;
	}
	SparseStructures::CSR save_csr;
	SparseStructures::CSR::readModFEMcrsMatrixFromFile(save_csr, load_file);
	load_file.close();
	size_t last_slash = file_name.rfind('/', file_name.size() - 1) + 1;
	std::string save_name = file_name.substr(last_slash, file_name.rfind(".txt") - last_slash) + "_binary";
	std::ofstream save_file(save_name, std::ifstream::binary);
	SparseStructures::CSR::saveModFEMcsrMatrixToBinary(save_csr, save_file);
	std::cout << "Done! Saved under:" << save_name << "\n";
	save_file.close();
#ifdef _DEBUG
	SparseStructures::CSR save2_csr;
	load_file.open(save_name, std::ifstream::binary);
	SparseStructures::CSR::readModFEMcsrMatrixFromBinary(save2_csr, load_file);
	//SparseStructures::CSR::print(save2_csr, 100);
	if (!(save2_csr == save_csr)) std::cout << "Houston, we have a problem.\n";
	std::cout << ((save2_csr == save_csr) ? "OK" : "NOT OK") << std::endl;
	load_file.close();
#endif
	return 0;
}

int main(int argc, char* argv[]) {
	
	std::string file_name;
	short file_type;
	short evaluation_type;
	short matrix_type;
	short load_ops;
	std::string solve_name;

	try {
		ArgumentParser parser(argc, argv);
#ifdef _DEBUG
		parser.print();
#endif
		file_name = parser.getArgument('f').first;
		file_type = parser.getArgument('l').second;
		evaluation_type = parser.getArgument('e').second;
		matrix_type = parser.getArgument('m').second;
		load_ops = parser.getArgument('t').second;
		solve_name = parser.getArgument('s').first;
	}
	catch (std::string error) {
		std::cout << error << std::endl;
		return -3;
	}
	
	std::cout << "Opening " << file_name << " " << std::endl;
	std::ifstream matrix_file(file_name);
	if (!matrix_file.is_open()) {
		std::cout << "Cannot open file " << file_name << "\n";
		return -1;
	}
	matrix_file.close();
	switch (load_ops) {
	case static_cast<short>(0):
		return  __cudssAndEigen(file_name, file_type, evaluation_type, matrix_type, solve_name);
	case static_cast<short>(1):
		std::cout << "Not implemented\n";
		return 0;
	case static_cast<short>(2):
		return __saveToBinary(file_name);
	case static_cast<short>(3):
		//{
		//	SparseStructures::CSR matrix;
		//	__loadMatrix(file_name, file_type, matrix);
		//	size_t perm_s = matrix.getN() * sizeof(int);
		//	//int* perm_col = nullptr;
		//	int* perm_col = new int[matrix.getN()];
		//	int* perm_row = new int[matrix.getN()+1];

		//	for (auto i = 0; i < matrix.getN(); i++) {
		//		perm_col[i] = perm_row[i] = i;
		//	}
		//	perm_row[matrix.getN()] = matrix.getN();
		//	saveSparcityImage("test.png", matrix.getN(), matrix.getNNZ(), matrix.getRowOff(), matrix.getColInd(), perm_col, perm_row, 2048);
		//	delete[] perm_col;
		//	delete[] perm_row;
		//}
	{
		int pk = 1024;
		if (!solve_name.empty()) pk = std::stoi(solve_name);
		return __getSparcityPatterns(file_name, file_type, matrix_type, pk);
	}
	default:
		break;
	}
	
	

	return 0;
}
