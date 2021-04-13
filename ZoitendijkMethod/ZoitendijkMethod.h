#pragma once
#include <vector>
#include "Linear.h"
#include "Simplex.h"

class ZoitendijkMethod {
	using xn_t = std::vector<double>;
public:
	ZoitendijkMethod(size_t dim, double lamda = 0.5);

	void set_limitaions(Limitations const& limitations);

	void calculate();

	~ZoitendijkMethod();

private:
	void init_first_approximation();

	static double dot_product(xn_t const& x, xn_t const& y);

	size_t dim;
	xn_t x0;
	Limitations limitations;
	double lamda;
};