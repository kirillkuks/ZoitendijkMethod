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
	double limitation_value(size_t index, xn_t const& x) const;
	xn_t solve_subtask(xn_t const& x, std::vector<size_t> const& almost_active);

	bool is_in_area(xn_t const& x) const;
	double find_next_alpha(xn_t const& xk, double eta_k, xn_t const& s_k);

	static double dot_product(xn_t const& x, xn_t const& y);
	static void scale(xn_t& vec, double mult);
	static xn_t add(xn_t const& x, xn_t const& y);

	std::vector<size_t> build_almost_active(xn_t const& x);
	double delta_null_active(xn_t const& x);

	size_t dim;
	xn_t x0;
	Function objective_function;
	Limitations limitations;
	std::vector<Function> gradients;
	double lamda;
	double delta;
};