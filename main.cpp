#include <iostream>
#include <fstream>
#include <random>

#include "cudss_wrapper.h"
#include "eigen_wrapper.h"
#include "matrixes/SparseStructures.h"
#include "ChronoTimer.h"

int main(int argc, char* argv[]) {

	const int max_args = 5;

	if (argc < 2) {
		std::cout << "Error: you need to provide file name.\n";
		return -4;
	}

	if (argc > max_args) {
		std::cout << "Correct usege: " << argv[0] << " "
			<< "<file_name> " 
			<< "<t for text file, b for binary> " 
			<< "<y - skip Eigen, n - don't skip> "
			<< "<y - skip Eigen, n - don't skip> "
			<< "<matrix type: 0 (General), 1 (Symetric), 2 (Hermitian)>"
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

	bool skip_eigen = false;
	if (argc >= 4 && argv[3][1] == '\0') {
		switch (argv[3][0]) {
		case 'y': case 'Y':
			skip_eigen = true;
			break;
		case 'n': case 'N':
			skip_eigen = false;
			break;
		default:
			std::cout << "Wrong option for comparison: y - skip, n - no skip\n";
			return -3;
		}
	}
	unsigned short matrix_type = 0;
	if (argc >= 5 && argv[4][1] == '\0') {
		switch (argv[4][0]) {
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
	std::string file_name(argv[1]);
	std::cout << "Opening " << file_name << " " << skip_eigen
		<< std::endl;

	std::ifstream matrix_file(file_name);
	if (!matrix_file.is_open()) {
		std::cout << "Cannot open file " << argv[1] << "\n";
		return -1;
	}

	SparseStructures::CSR test_csr;
	SparseStructures::CSR::readModFEMcrsMatrixFromFile(test_csr, matrix_file);

	matrix_file.close();

#ifdef _DEBUG
	SparseStructures::CSR::print(test_csr, 100);
#endif

	double* b = new double[test_csr.getN()];

	for (auto i = 0u; i < test_csr.getN(); i++)
		b[i] = 1.;

	double* x_1 = nullptr, * x_2 = nullptr;
	if (!skip_eigen)
		eigenDecompositon(test_csr, b, &x_2);

	size_t how_much = 10llu;
	ChronoTimer timer(6llu * how_much);
	for (auto i = 0llu; i < how_much; i++) {
		cuDSSDecompositionWithSynchronization(timer, test_csr, b, &x_1, 0, 0, 0);
		if (skip_eigen) {
			std::cout << "Skipped eigenDecomposition" << "\n";
			continue;
		}
		const size_t data_points = 7llu;
		double tolerance[data_points] = { 1.e-16, 1.e-15, 1.e-14, 1.e-13, 1.e-12, 1.e-11, 1.e-10 };
		int times[data_points] = { 0,0,0,0,0,0,0 };
		for (auto i = 0u; i < test_csr.getN(); i++)
			for (auto k = 0llu; k < data_points; k++)
				if (std::abs(x_1[i] - x_2[i]) < tolerance[k])
					times[k]++;

		std::cout << "For attempt " << i << "\n";
		double next_percent = static_cast<double>(times[0]) * 100. / static_cast<double>(test_csr.getN());
		size_t k = 0llu;
		do {
			std::cout << "Inside given tolerance: " << tolerance[k] << " " << times[k] << ", percent: " << next_percent << " % \n";
			if (next_percent >= 100.) break;
			next_percent = static_cast<double>(times[++k]) * 100. / static_cast<double>(test_csr.getN());
		} while (k != data_points);
	}

	size_t last_slash = file_name.rfind('/', file_name.size() - 1) + 1;
	std::string time_name = file_name.substr(last_slash, file_name.rfind(".txt") - last_slash) + "_test_time.txt";
	std::ofstream timefile(time_name);
	if (timefile.is_open()) {
		timer.saveTimesToFile(timefile, 6);
		timefile.close();
		std::cout << "Time stamps saved under file: " << time_name << "\n";
	}

	delete[] b;
	delete[] x_1;
	delete[] x_2;

	return 0;
}
