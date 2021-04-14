#pragma once
#include <vector>
#include <functional>

class Function {
	using xn_t = std::vector<double>;
public:
	Function(std::function<double(xn_t const&)> func, std::function<xn_t(xn_t const&)> grad);
	Function(std::vector<double> const& func, std::function<xn_t(xn_t const&)> grad);

	double operator()(xn_t const& x);
	xn_t gradient_at(xn_t const& x);

	static double dot_product(xn_t const& x, xn_t const& y);

	~Function();

private:
	enum class FuncType {
		NORMAL,
		LINEAR,
	};

	FuncType ft;
	std::function<double(xn_t const&)> func;
	std::vector<double> lin_func;
	std::function<xn_t(xn_t const&)> grad;

};
