#pragma once
#include <iostream>
#include <vector>
#include <limits>

#include "Matrix.h"

enum class KW {
	KW_LIMITATIONS = 0,
	KW_SIGN,
	KW_EQ,
	KW_GT,
	KW_LE,
	KW_NULL
};

enum class LT {
	LT_EQ = 0,
	LT_GT,
	LT_LE,
	LT_ERROR
};

struct Limitations {
public:
	Limitations();
	void add_limitations(std::pair<std::vector<double>, LT>&&);
	void add_limitations(std::pair<std::vector<double>, LT>&);
	size_t limititation_size() const;
	~Limitations();

	std::vector<std::pair<std::vector<double>, LT>> limitations;
};

struct Linear {
public:
	Linear(std::vector<double>&, Limitations&, std::vector<bool>&);

	//std::vector<double> solve_task();
	Matrix& get_matrix();
	std::vector<double>& get_b();
	std::vector<double>& get_obj_func();


	Linear* get_dual_program();
private:
	static size_t vars_in_canonical(std::vector<double>&, Limitations&, std::vector<bool>&);
	//static bool next_combination(std::vector<size_t>&, size_t);

	Linear(std::vector<double>&, Limitations&, std::vector<bool>&, Linear*);
	void to_canonical(std::vector<double>&, Limitations&, std::vector<bool>&);
	void to_dual(std::vector<double>&, Limitations&, std::vector<bool>&);

	//bool is_linear_independence(Matrix&);
	double calculate_objective(std::vector<double>&);
	//Matrix sub_matrix(std::vector<size_t>&);
	//double determinant(Matrix);
	//bool in_allowable_area(std::vector<double>&);
	//std::vector<double> back_to_original_vars(std::vector<double>&);

	size_t original_dimension;

	Matrix A;
	std::vector<double> b;
	std::vector<double> objective_function;
	TT task_type;

	std::vector<std::pair<size_t, size_t>> original_vars;

	Linear* dual_program;
};
