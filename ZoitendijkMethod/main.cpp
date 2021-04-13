#include <iostream>
#include "Linear.h"
#include "ZoitendijkMethod.h"

int main() {
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
	limitations.add_limitations({ {2, 1, 2}, LT::LT_LE });
	limitations.add_limitations({ {-1, 1, 3}, LT::LT_LE });
	limitations.add_limitations({ {3, -2, 4}, LT::LT_LE });
	limitations.add_limitations({ {1, 1, -5}, LT::LT_GT });
	ZoitendijkMethod method(2);
	method.set_limitaions(limitations);
	method.calculate();
	return 0;
}