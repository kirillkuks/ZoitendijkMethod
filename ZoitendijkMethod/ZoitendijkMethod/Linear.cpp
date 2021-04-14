#include "Linear.h"

Limitations::Limitations() {}

void Limitations::add_limitations(std::pair<std::vector<double>, LT>&& limitation) {
	limitations.push_back(limitation);
}

void Limitations::add_limitations(std::pair<std::vector<double>, LT>& limitation) {
	limitations.push_back(limitation);
}

size_t Limitations::limititation_size() const {
	return limitations.size();
}

Limitations::~Limitations() {}


Linear::Linear(std::vector<double>& function, Limitations& limitations, std::vector<bool>& vars_sign)
	: A(limitations.limitations.size(), vars_in_canonical(function, limitations, vars_sign)),
	b(limitations.limitations.size()),
	objective_function(A.get_m()),
	original_dimension(function.size()),
	original_vars(function.size()),
	task_type{ TT::TT_MIN },
	dual_program(nullptr) {

	to_canonical(function, limitations, vars_sign);

	to_dual(function, limitations, vars_sign);

	/*std::cout << "\n*****************************************\n";
	for (auto elem : function) {
		std::cout << elem << " ";
	}
	std::cout << "\n#########################################\n";
	for (size_t i = 0; i < limitations.limitations.size(); ++i) {
		for (size_t j = 0; j < limitations.limitations[i].first.size(); ++j) {
			std::cout << limitations.limitations[i].first[j] << " ";
		}
		std::cout << " // sign: " << limitations.limitations[i].second << std::endl;
	}
	std::cout << "#########################################\n";
	for (auto elem : vars_sign) {
		std::cout << elem << " ";
	}
	std::cout << "\n*****************************************\n";*/

	dual_program = new Linear(function, limitations, vars_sign, this);
}

Linear::Linear(std::vector<double>& function, Limitations& limitations, std::vector<bool>& vars_sign, Linear* dual_program)
	: A(limitations.limitations.size(), vars_in_canonical(function, limitations, vars_sign)),
	b(limitations.limitations.size()),
	objective_function(A.get_m()),
	original_dimension(function.size()),
	original_vars(function.size()),
	task_type{ TT::TT_MAX },
	dual_program{ dual_program } {

	to_canonical(function, limitations, vars_sign);
}

void Linear::to_dual(std::vector<double>& function, Limitations& limitations, std::vector<bool>& vars_sign) {
	std::vector<double> dual_function(limitations.limitations.size());
	Limitations dual_limitations;
	std::vector<bool> dual_vars_sign(dual_function.size());

	for (size_t i = 0; i < dual_function.size(); ++i) {
		if (limitations.limitations[i].second == LT::LT_LE) {
			for (size_t j = 0; j < original_dimension + 1; ++j) {
				limitations.limitations[i].first[j] *= -1;
			}
		}
	}

	for (size_t i = 0; i < dual_function.size(); ++i) {
		dual_function[i] = limitations.limitations[i].first[original_dimension];
	}

	for (size_t i = 0; i < original_dimension; ++i) {
		std::vector<double> new_limitation(dual_function.size() + 1);
		for (size_t j = 0; j < dual_function.size(); ++j) {
			new_limitation[j] = limitations.limitations[j].first[i];
		}
		new_limitation[dual_function.size()] = function[i];
		dual_limitations.add_limitations({ new_limitation, vars_sign[i] ? LT::LT_LE : LT::LT_EQ });
	}

	for (size_t i = 0; i < dual_function.size(); ++i) {
		if (limitations.limitations[i].second != LT::LT_EQ) {
			dual_vars_sign[i] = true;
		}
	}

	std::swap(dual_function, function);
	std::swap(dual_limitations, limitations);
	std::swap(dual_vars_sign, vars_sign);
}

void Linear::to_canonical(std::vector<double>& function, Limitations& limitations, std::vector<bool>& vars_sign) {
	for (size_t i = 0; i < b.size(); ++i) {
		b[i] = limitations.limitations[i].first[function.size()];
	}

	size_t m = 0;
	for (size_t i = 0; i < function.size(); ++i, ++m) {
		if (vars_sign[i]) {
			objective_function[m] = function[i];
			original_vars[i] = { m, m };
		}
		else {
			objective_function[m++] = function[i];
			objective_function[m] = -function[i];
			original_vars[i] = { m - 1, m };
		}
	}

	for (size_t i = 0, n = m; i < limitations.limitations.size(); ++i) {

		for (size_t k = 0, j = 0; k < limitations.limitations[i].first.size() - 1; ++k, ++j) {
			if (vars_sign[k]) {
				A[i][j] = limitations.limitations[i].first[k];
			}
			else {
				A[i][j++] = limitations.limitations[i].first[k];
				A[i][j] = -limitations.limitations[i].first[k];
			}
		}

		if (limitations.limitations[i].second == LT::LT_GT) {
			A[i][n++] = -1;
		}
		if (limitations.limitations[i].second == LT::LT_LE) {
			A[i][n++] = 1;
		}

	}

	//Для проверки
	std::cout << "c:\n";
	for (auto elem : objective_function) {
		std::cout << elem << ' ';
	}
	std::cout << "\nA:\n";
	A.print();
	std::cout << "b:\n";
	for (auto elem : b) {
		std::cout << elem << ' ';
	}
	std::cout << std::endl;
}

Linear* Linear::get_dual_program() {
	return dual_program;
}

size_t Linear::vars_in_canonical(std::vector<double>& function, Limitations& limitations, std::vector<bool>& vars_sign) {
	size_t vars_num = function.size();
	for (auto elem : vars_sign) {
		if (!elem) {
			++vars_num;
		}
	}
	for (auto vec : limitations.limitations) {
		if (vec.second != LT::LT_EQ) {
			++vars_num;
		}
	}
	return vars_num;
}

double Linear::calculate_objective(std::vector<double>& x) {
	double res = 0;
	for (size_t i = 0; i < objective_function.size(); ++i) {
		res += objective_function[i] * x[i];
	}
	return res;
}

Matrix& Linear::get_matrix() {
	for (size_t i = 0; i < A.get_n(); ++i) {
		for (size_t j = 0; j < A.get_m(); ++j) {
			std::cout << A[i][j] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << "End\n";
	return A;
}

std::vector<double>& Linear::get_b() {
	return b;
}

std::vector<double>& Linear::get_obj_func() {
	return objective_function;
}

/*bool Linear::next_combination(std::vector<size_t>& vec, size_t n) {
	size_t k = vec.size();
	int i;
	for (i = k - 1; i >= 0 && vec[i] == n - 1; --i, --n) {}
	if (i < 0) {
		return false;
	}
	++vec[i];
	for (size_t j = i + 1; j < k; ++j) {
		vec[j] = vec[j - 1] + 1;
	}
	return true;
}

double Linear::determinant(Matrix matrix) {
	size_t size = b.size();
	double det = 1;
	for (size_t i = 0; i < size; ++i) {
		size_t k = i;
		for (size_t j = i + 1; j < size; ++j) {
			if (std::abs(matrix[j][i]) > std::abs(matrix[k][i])) {
				k = j;
			}
		}
		if (std::abs(matrix[k][i]) < 1E-5) {
			return 0;
		}
		for (size_t j = 0; j < size; ++j) {
			std::swap(matrix[i][j], matrix[k][j]);
		}
		if (i != k) {
			det = -det;
		}
		det *= matrix[i][i];
		for (size_t j = i + 1; j < size; ++j) {
			matrix[i][j] /= matrix[i][i];
		}
		for (size_t j = 0; j < size; ++j) {
			if (j != i && std::abs(matrix[j][i]) > 1E-5) {
				for (size_t k = i + 1; k < size; ++k) {
					matrix[j][k] -= matrix[i][k] * matrix[j][i];
				}
			}
		}
	}
	return det;
}

bool Linear::is_linear_independence(Matrix& submatrix) {
	return determinant(submatrix) != 0;
}

bool Linear::in_allowable_area(std::vector<double>& vec) {
	for (auto elem : vec) {
		if (elem < 0) {
			return false;
		}
	}
	return true;
}

Matrix Linear::sub_matrix(std::vector<size_t>& vec) {
	size_t size = vec.size();
	Matrix submatrix(size, size);
	for (size_t j = 0; j < size; ++j) {
		for (size_t i = 0; i < size; ++i) {
			submatrix[i][j] = A[i][vec[j]];
		}
	}
	return submatrix;
}

std::vector<double> Linear::back_to_original_vars(std::vector<double>& x) {
	std::cout << "\nCanonical: ";
	for (auto elem : x) {
		std::cout << elem << ' ';
	}
	std::cout << std::endl;
	std::vector<double> vars(original_dimension);
	for (size_t i = 0; i < original_dimension; ++i) {
		if (original_vars[i].first == original_vars[i].second) {
			vars[i] = x[original_vars[i].first];
		}
		else {
			vars[i] = x[original_vars[i].first] - x[original_vars[i].second];
		}
	}
	return vars;
}

std::vector<double> Linear::solve_task() {
	size_t canonical_dimension = objective_function.size();
	std::vector<double> optimal(canonical_dimension);
	double limit_value = task_type == TT::TT_MIN ? std::numeric_limits<double>::infinity() : -std::numeric_limits<double>::infinity();

	std::vector<size_t> vectors_in_basis(b.size());
	for (size_t i = 0; i < vectors_in_basis.size(); ++i) {
		vectors_in_basis[i] = i;
	}

	do {
		Matrix submatrix = sub_matrix(vectors_in_basis);
		if (is_linear_independence(submatrix)) {
			std::vector<double> x(canonical_dimension);
			std::vector<double> res(b.size());

			res = SeidelMethod::solve(submatrix, b);

			if (in_allowable_area(res)) {
				for (size_t i = 0; i < b.size(); ++i) {
					x[vectors_in_basis[i]] = res[i];
				}

				std::cout << "x:\n";
				for (auto elem : x) {
					std::cout << elem << " ";
				}
				std::cout << std::endl;

				double potentional_optimal = calculate_objective(x);
				if ((task_type == TT::TT_MIN && limit_value > potentional_optimal) || (task_type == TT::TT_MAX && limit_value < potentional_optimal)) {
					optimal = x;
					limit_value = potentional_optimal;
				}

			}
		}
	} while (next_combination(vectors_in_basis, canonical_dimension));

	std::cout << "\nMin: " << min;
	std::cout << "\nVector:\n";
	for (auto elem : optimal) {
		std::cout << elem << ' ';
	}
	std::vector<double> vars = back_to_original_vars(optimal);
	std::cout << "Original:\n";
	for (auto elem : vars) {
		std::cout << elem << ' ';
	}

	return back_to_original_vars(optimal);
}*/