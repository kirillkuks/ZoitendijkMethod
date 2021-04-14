#include "SeidelMethod.h"

std::vector<double> SeidelMethod::solve(Matrix& A, std::vector<double>& b) {
	size_t size = b.size();
	std::vector<double> x(size), g(size), x0(size), err(size);
	Matrix B = A.get_transposed();
	Matrix C = multy(B, A);
	g = b_create(B, b);
	size_t counter = 0;
	do {
		for (size_t j = 0; j < x0.size(); ++j) {
			x0[j] = x[j];
		}
		for (size_t i = 0; i < size; ++i) {
			double sum = 0;
			for (size_t j = 0; j < i; ++j) {
				sum += (C[i][j] * x[j]);
			}
			for (size_t j = i + 1; j < size; ++j) {
				sum += (C[i][j] * x0[j]);
			}
			x[i] = (g[i] - sum) / C[i][i];
		}
		err = b_create(A, x);
		for (size_t j = 0; j < size; ++j) {
			err[j] -= b[j];
		}
		//std::cout << "Norm: " << norm(err) << std::endl;
		++counter;
	} while (norm(err) > 1E-5 && counter < 200);
	return x;
}

Matrix SeidelMethod::multy(Matrix& A, Matrix& B) {
	Matrix C(A.get_n(), B.get_m());
	for (size_t i = 0; i < A.get_n(); ++i) {
		for (size_t j = 0; j < B.get_m(); ++j) {
			C[i][j] = 0;
			for (size_t k = 0; k < A.get_m(); ++k) {
				C[i][j] += A[i][k] * B[k][j];
			}
		}
	}
	return C;
}

std::vector<double> SeidelMethod::b_create(Matrix& A, std::vector<double>& x) {
	std::vector<double> b(x.size());
	for (size_t i = 0; i < b.size(); ++i) {
		b[i] = 0;
		for (size_t j = 0; j < b.size(); ++j) {
			b[i] += A[i][j] * x[j];
		}
	}
	return b;
}

double SeidelMethod::norm(std::vector<double>& x) {
	double res = 0;
	for (auto elem : x) {
		res += elem * elem;
	}
	return std::sqrt(res);
}