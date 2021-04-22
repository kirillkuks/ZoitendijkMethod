#include <iostream>
#include <cassert>
#include "Linear.h"
#include "ZoitendijkMethod.h"

#define N 2

using xn_t = std::vector<double>;

int main() {
	auto a = 2, b = 4, c = 5;
	auto func = [a, b, c](xn_t const& x) -> double {
		assert(x.size() == N);
		auto x1 = x[0], x2 = x[1];
		return a * x1 + x2 + 4 * sqrt(1 + b * x1 * x1 + c * x2 * x2);
	};
	auto grad = [a, b, c](xn_t const& x) -> xn_t {
		auto size = x.size();
		assert(size == N);

		auto x1 = x[0], x2 = x[1];
		auto value = sqrt(1 + b * x1 * x1 + c * x2 * x2);
		return { a + 4.0 * b * x1 / value, 1 + 4.0 * c * x2 / value };
	};

	/*std::vector<double> of = { 1, 0, 0 };
	Limitations limitaions;
	limitaions.add_limitations({ {0, 1, 0, -1}, LT::LT_GT });
	limitaions.add_limitations({ {0, 0, 1, -1}, LT::LT_GT });
	limitaions.add_limitations({ {0, 1, 0, 1}, LT::LT_LE });
	limitaions.add_limitations({ {0, 0, 1, 1}, LT::LT_LE });
	limitaions.add_limitations({ {-1, 0, -8, 0}, LT::LT_LE });
	std::vector<bool> signs = { false, false, false };
	Linear linear(of, limitaions, signs);*/

	Limitations limitations;

	// Точка внутри области
	/*limitations.add_limitations({ {2, 1, 2}, LT::LT_LE });
	std::function<xn_t(xn_t const&)> grad1 = [](xn_t const&) -> xn_t { return { 2, 1 }; };
	Function func1([](xn_t const& x) -> double { return 2 * x[0] + x[1] - 2; }, grad1);

 	limitations.add_limitations({ {-1, 1, 3}, LT::LT_LE });
	std::function<xn_t(xn_t const&)> grad2 = [](xn_t const&) -> xn_t { return { -1, 1 }; };
	Function func2([](xn_t const& x) -> double { return - x[0] + x[1] - 3; }, grad2);

	limitations.add_limitations({ {3, -2, 4}, LT::LT_LE });
	std::function<xn_t(xn_t const&)> grad3 = [](xn_t const&) -> xn_t { return { 3, -2 }; };
	Function func3([](xn_t const& x) -> double { return 3 * x[0] - 2 * x[1] - 4; }, grad3);

	limitations.add_limitations({ {-1, -1, 5}, LT::LT_LE });
	std::function<xn_t(xn_t const&)> grad4 = [](xn_t const&) -> xn_t { return { -1, -1 }; };
	Function func4([](xn_t const& x) -> double { return - x[0] - x[1] - 5; }, grad4);*/

	// Точка на границе
	limitations.add_limitations({ {2, 1, 2}, LT::LT_LE });
	std::function<xn_t(xn_t const&)> grad1 = [](xn_t const&) -> xn_t { return { 2, 1 }; };
	Function func1([](xn_t const& x) -> double { return 2 * x[0] + x[1] -2; }, grad1);

	limitations.add_limitations({ {-1, 1, 3}, LT::LT_LE });
	std::function<xn_t(xn_t const&)> grad2 = [](xn_t const&) -> xn_t { return { -1, 1 }; };
	Function func2([](xn_t const& x) -> double { return -x[0] + x[1] - 3; }, grad2);

	limitations.add_limitations({ {3, -2, 4}, LT::LT_LE });
	std::function<xn_t(xn_t const&)> grad3 = [](xn_t const&) -> xn_t { return { 3, -2 }; };
	Function func3([](xn_t const& x) -> double { return 3 * x[0] - 2 * x[1] -4; }, grad3);

	limitations.add_limitations({ {-1, -1, 0.182}, LT::LT_LE });
	std::function<xn_t(xn_t const&)> grad4 = [](xn_t const&) -> xn_t { return { -1, -1 }; };
	Function func4([](xn_t const& x) -> double { return -x[0] - x[1] - 0.182; }, grad4);

	Function function(func, grad);
	ZoitendijkMethod method(function, 2);

	method.set_limitaions(limitations, { func1, func2, func3, func4 });
	method.calculate();
	return 0;
}