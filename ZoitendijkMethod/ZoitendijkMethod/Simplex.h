#pragma once

#include <vector>
#include <iostream>
#include "Matrix.h"


class Simplex {
public:
	Simplex(Matrix&, std::vector<double>&, std::vector<double>&, TT type_task);
	std::vector<double> answer_func();

private:
	std::vector<std::pair<std::vector<double>, size_t>> data;
	std::vector<double> func;
	std::vector<double> delta;
	double answer;
	bool have_ans;

	void Choose(size_t i, size_t min);
	size_t Determine(size_t i);
	size_t Check();
	size_t Num_Var(std::vector<double>&);
	bool Check_Data();
	void Positive_b();
};
