#pragma once
#include <vector>

class Solver
{
public:
	static std::vector<double> solveGaussSimple(const std::vector<std::vector<double>>& matrix, const std::vector<double>& results);
};

