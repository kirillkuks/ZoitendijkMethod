#include "Matrix.h"

Matrix::Matrix(size_t n, size_t m = 0) : n{ n }, m{ m == 0 ? n : m } {
	matrix = new double* [n];
	for (size_t i = 0; i < n; ++i) {
		matrix[i] = new double[m];
		for (size_t j = 0; j < m; ++j) {
			matrix[i][j] = 0;
		}
	}
}

Matrix::Matrix(Matrix const& other) : Matrix(other.n, other.m) {
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = 0; j < m; ++j) {
			matrix[i][j] = other.matrix[i][j];
		}
	}
}

Matrix::Vector::Vector(double* vector) : vector{ vector } {}

Matrix& Matrix::operator=(Matrix other) {
	std::swap(matrix, other.matrix);
	std::swap(n, other.n);
	std::swap(m, other.m);
	return *this;
}

double& Matrix::Vector::operator[](size_t index) {
	return vector[index];
}

Matrix::Vector Matrix::operator[](size_t index) {
	return Vector(matrix[index]);
}

std::vector<double> Matrix::gauss(std::vector<double>& b) {
	for (size_t i = 0; i < n; i++) {
		size_t str = 0;
		for (size_t j = 0; j < m; j++) {
			if (fabs(matrix[i][j]) > 0.00001) {
				str = j;
				break;
			}
		}
		for (size_t j = 0; j < m; j++) {
			if (j != str)
				matrix[i][j] /= matrix[i][str];
		}
		b[i] /= matrix[i][str];
		matrix[i][str] = 1.;
		for (size_t j = 0; j < n; j++) {
			if (j != i) {
				for (size_t q = 0; q < m; q++) {
					if (q != str)
						matrix[j][q] -= matrix[j][str] * matrix[i][q];
				}
				b[j] -= matrix[j][str] * b[i];
				matrix[j][str] = 0.;
			}
		}
	}
	return b;
}

void Matrix::Change_Str(size_t i, size_t j) {
	for (size_t q = 0; q < m; q++) {
		double num = matrix[i][q];
		matrix[i][q] = matrix[j][q];
		matrix[j][q] = num;
	}
}

void Matrix::Change_Col(size_t i, size_t j) {
	for (size_t q = 0; q < n; q++) {
		double num = matrix[q][i];
		matrix[q][i] = matrix[q][j];
		matrix[q][j] = num;
	}
}

void Matrix::print() const {
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = 0; j < m; ++j) {
			std::cout << matrix[i][j] << ' ';
		}
		std::cout << std::endl;
	}
}

size_t Matrix::get_n() const {
	return n;
}

size_t Matrix::get_m() const {
	return m;
}

Matrix Matrix::get_transposed() const {
	Matrix AT(m, n);
	for (size_t i = 0; i < m; ++i) {
		for (size_t j = 0; j < n; ++j) {
			AT[i][j] = matrix[j][i];
		}
	}
	return AT;
}

Matrix::~Matrix() {
	for (size_t i = 0; i < n; ++i) {
		delete[] matrix[i];
	}
	delete[] matrix;
}