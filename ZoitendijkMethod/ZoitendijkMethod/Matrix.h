#pragma once

#include <iostream>
#include <vector>

enum TT {
	TT_MIN = 0,
	TT_MAX
};

struct Matrix {
private:
	struct Vector;
public:
	Matrix(size_t, size_t);
	Matrix(Matrix const&);

	Matrix& operator=(Matrix);
	Vector operator[](size_t);

	size_t get_n() const;
	size_t get_m() const;
	Matrix get_transposed() const;
	std::vector<double> gauss(std::vector<double>&);
	void Change_Str(size_t i, size_t j);
	void Change_Col(size_t i, size_t j);

	void print() const;

	~Matrix();
private:
	struct Vector {
		Vector(double*);
		double& operator[](size_t);
	private:
		double* vector;
	};

	double** matrix;
	size_t n, m;
};
