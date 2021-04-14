#pragma once
#include <vector>
#include "Linear.h"
#include "Simplex.h"
#include "Function.h"

class ZoitendijkMethod {
	using xn_t = std::vector<double>;
public:
	ZoitendijkMethod(Function const& function, size_t dim, double lamda = 0.5);

	void set_limitaions(Limitations const& limitations, std::vector<Function> const& grads);
	void set_function(Function const& objective_function);

	void calculate();

	~ZoitendijkMethod();

private:
	void init_first_approximation();
	double limitation_value(size_t index, xn_t const& x);

	static double dot_product(xn_t const& x, xn_t const& y);
	std::vector<size_t> build_almost_active(xn_t const& x);

	size_t dim;
	xn_t x0;
	Function objective_function;
	Limitations limitations;
	std::vector<Function> gradients;
	double lamda;
	double delta;
};