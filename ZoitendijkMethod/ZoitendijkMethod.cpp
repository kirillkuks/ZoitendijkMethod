#include <cassert>

#include "ZoitendijkMethod.h"

double ZoitendijkMethod::dot_product(xn_t const& x, xn_t const& y) {
	size_t size = x.size();
	assert(size == y.size());
	
	double val = 0;
	for (size_t i = 0; i < size; ++i) {
		val += x[i] * y[i];
	}
	return val;
}

ZoitendijkMethod::ZoitendijkMethod(size_t dim, double lamda) : dim{ dim }, limitations {  }, x0{  }, lamda{ lamda } {}

ZoitendijkMethod::~ZoitendijkMethod() {}

void ZoitendijkMethod::set_limitaions(Limitations const& limitations_) {
	limitations = limitations_;
	for (auto& limitation : limitations.limitations) {
		if (limitation.second == LT::LT_GT) {
			for (auto& elem : limitation.first) {
				elem *= -1;
			}
			limitation.second = LT::LT_LE;
		}
	}
}

void ZoitendijkMethod::init_first_approximation() {
	std::vector<double> objective_function(dim + 1);
	objective_function[0] = 1;
	Limitations lims;
	for (auto const& limitation : limitations.limitations) {
		std::pair<std::vector<double>, LT> cond = std::make_pair(std::vector<double>(dim + 2), LT::LT_LE);
		cond.first[0] = -1;
		for (size_t i = 0; i < dim + 1; ++i) {
			cond.first[i + 1] = limitation.first[i];
		}
		lims.add_limitations(cond);
	}

	for (auto& lim : lims.limitations) {
		for (auto elem : lim.first) {
			std::cout << elem << " | ";
		}
		std::cout << std::endl;
	}
	std::vector<bool> var_signs(dim + 1);
	Linear subtask(objective_function, lims, var_signs);
	Simplex simplex(subtask.get_matrix(), subtask.get_b(), subtask.get_obj_func(), TT::TT_MIN);
	auto optimal = simplex.answer_func();
	optimal = subtask.back_to_original_vars(optimal);
	std::cout << "Optimal" << std::endl;
	for (size_t i = 0; i < optimal.size(); i++) {
		std::cout << optimal[i] << " ";
	}
	std::cout << std::endl;
}

void ZoitendijkMethod::calculate() {
	init_first_approximation();
}
