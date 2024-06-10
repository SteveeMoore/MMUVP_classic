#include "Solver.h"
#include <stdexcept>

std::vector<double> Solver::solveGaussSimple(const std::vector<std::vector<double>>& matrix, const std::vector<double>& results)
{
    int n = results.size();
    std::vector<std::vector<double>> mat = matrix;
    std::vector<double> res = results;
    std::vector<double> x(n, 0.0);

    // Прямой ход метода Гаусса
    for (int i = 0; i < n; i++) {
        if (mat[i][i] == 0.0) {
            throw std::runtime_error("Нулевой элемент на диагонали, нельзя решить методом Гаусса без перестановок.");
        }

        for (int k = i + 1; k < n; k++) {
            double c = -mat[k][i] / mat[i][i];
            for (int j = i; j < n; j++) {
                mat[k][j] += c * mat[i][j];
            }
            res[k] += c * res[i];
        }
    }

    // Обратный ход метода Гаусса
    for (int i = n - 1; i >= 0; i--) {
        x[i] = res[i] / mat[i][i];
        for (int k = i - 1; k >= 0; k--) {
            res[k] -= mat[k][i] * x[i];
        }
    }
    return x;
}
