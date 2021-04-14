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

ZoitendijkMethod::ZoitendijkMethod(Function const& function, size_t dim, double lamda) 
	: objective_function{function}, dim{ dim }, limitations{  }, x0{ std::vector<double>(dim) }, lamda{ lamda }, delta{ 0 } {}

ZoitendijkMethod::~ZoitendijkMethod() {}

void ZoitendijkMethod::set_limitaions(Limitations const& limitations_, std::vector<Function> const& gradients_) {
	//std::cout << limitations.limitations.size() << " <> " << gradients_.size();
	assert(limitations_.limititation_size() == gradients_.size());
	limitations = limitations_;
	gradients = gradients_;
	for (auto& limitation : limitations.limitations) {
		if (limitation.second == LT::LT_GT) {
			for (auto& elem : limitation.first) {
				elem *= -1;
			}
			limitation.second = LT::LT_LE;
		}
	}
}

void ZoitendijkMethod::set_function(Function const& objective_func_) {
	objective_function = objective_func_;
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
	std::vector<double> optimal = simplex.answer_func();
	optimal = subtask.back_to_original_vars(optimal);
	//std::vector<double> optimal = subtask.solve_task();
	delta = -optimal[0];
	for (size_t i = 0; i < dim; ++i) {
		x0[i] = optimal[i + 1];
	}
	std::cout << "Delta: " << delta << std::endl;
	std::cout << "First Approximation\n";
	for (auto elem : x0) {
		std::cout << elem << ", ";
	}
	std::cout << std::endl;
	//Simplex simplex(subtask.get_matrix(), subtask.get_b(), subtask.get_obj_func(), TT::TT_MIN);
	//simplex.answer_func();
}

void ZoitendijkMethod::calculate() {
	init_first_approximation();
	for (size_t i = 0; i < limitations.limititation_size(); ++i) {
		std::cout << "Value" << i << ": " << limitation_value(i, x0) << std::endl;
	}
	auto almost_active = build_almost_active(x0);

	{
		Limitations new_limitations;
		/*auto limit = std::make_pair(std::vector<double>(dim + 2), LT::LT_LE);
		limit.first[0] = -1;
		xn_t grad = objective_function.gradient_at(x0);
		for (size_t i = 0; i < dim; ++i) {
			limit.first[i + 1] = grad[i];
		}
		limit.first[dim + 1] = 0;
		new_limitations.add_limitations(limit);*/
		for (size_t i = 0; i < almost_active.size(); ++i) {
			std::pair<std::vector<double>, LT> lim = std::make_pair(std::vector<double>(dim + 2), LT::LT_LE);
			lim.first[0] = -1;
			xn_t grad_val = gradients[almost_active[i]].gradient_at(x0);
			for (size_t j = 0; j < dim; ++j) {
				lim.first[j + 1] = grad_val[j];
			}
			lim.first[dim + 1] = 0;
			new_limitations.add_limitations(lim);
		}
		for (size_t i = 0; i < dim; ++i) {
			std::pair<std::vector<double>, LT> lim1 = std::make_pair(std::vector<double>(dim + 2), LT::LT_LE);
			std::pair<std::vector<double>, LT> lim2 = std::make_pair(std::vector<double>(dim + 2), LT::LT_GT);
			lim1.first[i + 1] = lim2.first[i + 1] = 1;
			lim1.first[dim + 1] = 1;
			lim2.first[dim + 1] = -1;
			new_limitations.add_limitations(lim1);
			new_limitations.add_limitations(lim2);
		}
		std::cout << "New Lims:\n";
		for (auto& lim : new_limitations.limitations) {
			for (auto elem : lim.first) {
				std::cout << elem << " | ";
			}
			std::cout << (lim.second == LT::LT_LE ? "LE" : "GT") << std::endl;
		}
		std::vector<double> obj_func(dim + 1);
		obj_func[0] = 1;
		std::vector<bool> var_signs(dim + 1);
		Linear linear(obj_func, new_limitations, var_signs);
		Simplex simplex(linear.get_matrix(), linear.get_b(), linear.get_obj_func(), TT::TT_MIN);
		std::vector<double> optimal = simplex.answer_func();
		std::cout << "OPT:\n";
		for(auto elem : optimal) {
			std::cout << elem << " ";
		}
		optimal = linear.back_to_original_vars(optimal);
		//linear.solve_task();

		std::cout << "Optimal:\n";
		for (auto elem : optimal) {
			std::cout << elem << " ";
		}
		std::cout << std::endl;
	}

}

std::vector<size_t> ZoitendijkMethod::build_almost_active(xn_t const& x) {
	std::vector<size_t> almost_active;
	for (size_t i = 0; i < limitations.limititation_size(); ++i) {
		double val = limitation_value(i, x);
		if (-delta <= val && val <= 0) {
			almost_active.push_back(i);
			std::cout << "Number: " << i << std::endl;
		}
	}
	return almost_active;
}

double ZoitendijkMethod::limitation_value(size_t index, xn_t const& x) {
	size_t size = x.size();
	assert(index < limitations.limititation_size());

	double val = 0;
	for (size_t i = 0; i < size; ++i) {
		val += limitations.limitations[index].first[i] * x[i];
	}
	val += -limitations.limitations[index].first[size];
	return val;
}
