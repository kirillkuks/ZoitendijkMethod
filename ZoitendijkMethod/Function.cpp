#include <cassert>

#include "Function.h"

Function::Function(std::function<double(xn_t const&)> func, std::function<xn_t(xn_t const&)> grad) : func{ func }, grad{ grad }, ft{ FuncType::NORMAL } {}

Function::Function(std::vector<double> const& func, std::function<xn_t(xn_t const&)> grad) : lin_func{ func }, grad{ grad }, ft{ FuncType::LINEAR } {}

Function::~Function() {}

double Function::dot_product(xn_t const& x, xn_t const& y) {
	size_t size = x.size();
	assert(size == y.size());

	double val = 0;
	for (size_t i = 0; i < size; ++i) {
		val += x[i] * y[i];
	}
	return val;
}

double Function::operator()(xn_t const& x) {
	switch (ft) {
	case FuncType::NORMAL:
		return func(x);
	case FuncType::LINEAR:
		return dot_product(lin_func, x);
	}
	return 0;
}

std::vector<double> Function::gradient_at(xn_t const& x) {
	return grad(x);
}