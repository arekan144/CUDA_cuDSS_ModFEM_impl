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

static int __getSparcityPatterns(std::string& file_name, short file_type, short matrix_type) {
	SparseStructures::CSR test_csr;
	__loadMatrix(file_name, file_type, test_csr);
	double* b = new double[test_csr.getN()];

	for (auto i = 0; i < test_csr.getN(); i++)
		b[i] = 1.;
	double* x_1 = nullptr, * x_2 = nullptr;
	cuDSSOnlyAnalisysAndSpPattern(test_csr, b, &x_1, matrix_type, 0, 0);
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
	

	/*if (!file_name_set) {
		std::cout << "Error: you need to provide file name.\n";
		return -256;
	}*/

	/*
	constexpr int max_args = 6;

	if (argc < 2) {
		std::cout << "Error: you need to provide file name.\n";
		return -4;
	}
	if (argc == 3 && argv[1][1] == '\0' && argv[1][0] == 's')
	{
		std::string file_name(argv[2]);
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
		if (!(save2_csr == save_csr)) std::cout << "Houston, we have a problem\n";
		std::cout << ((save2_csr == save_csr)?"OK":"NOT OK") << std::endl;
		load_file.close();
#endif
		return 0;
	}
	if (argc == 3 && argv[1][1] == '\0' && argv[1][0] == '')
	{

	}
	else if (argc == 3 && argv[1][1] == '\0' && argv[1][0] != 's')
	{
		std::cout << "Save feauture expects 's' as first argument, then 'file_name'\n";
		return -2;
	}

	if (argc > max_args) {
		std::cout << "Correct usege: " << argv[0] << " "
			<< "<file_name> " 
			<< "<t for text file, b for binary, OPT def= t> " 
			<< "<matrix type: 0 (General), 1 (Symetric), 2 (Hermitian), OPT def= 0> "
			<< "<y - skip Eigen, s - save the resoults, n - don't skip, l - load solution from file (adv skip), OPT def= n> "
			<< "<file_name_with_prev_result, OPT, ignored if prev != l> "
			<< "\n";
		return -2;
	}

	char load_ops = 't';
	if (argc >= 3 && argv[2][1] == '\0') {
		switch (argv[2][0]) {
		case 't': case 'T':
			load_ops = 't';
			break;
		case 'b': case 'B':
			load_ops = 'b';
			break;
		default:
			std::cout << "Wrong option: t for text file, b for binary\n";
			return -3;
		}
	}
	unsigned short matrix_type = 0;
	if (argc >= 4 && argv[3][1] == '\0') {
		switch (argv[3][0]) {
		case '0':
			matrix_type = static_cast<unsigned short>(0u);
			break;
		case '1':
			matrix_type = static_cast<unsigned short>(1u);
			break;
		case '2':
			matrix_type = static_cast<unsigned short>(2u);
			break;
		default:
			std::cout << "Wrong option for matrix type: 0, 1, 2\n";
			return -3;
		}
	}
	unsigned short nskip_eigen = static_cast<unsigned short>(1u);
	if (argc >= 5 && argv[4][1] == '\0') {
		switch (argv[4][0]) {
		case 'y': case 'Y':
			nskip_eigen = static_cast<unsigned short>(0u);
			break;
		case 'n': case 'N':
			nskip_eigen = static_cast<unsigned short>(1u);
			break;
		case 's': case 'S':
			nskip_eigen = static_cast<unsigned short>(2u);
			break;
		case 'l': case 'L':
			nskip_eigen = static_cast<unsigned short>(3u);
			break;
		default:
			std::cout << "Wrong option for comparison: y - skip, n - no skip, s - save\n";
			return -3;
		}
	}
	*/
	
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
		return __getSparcityPatterns(file_name, file_type, matrix_type);
	default:
		break;
	}
	
	

	return 0;
}
