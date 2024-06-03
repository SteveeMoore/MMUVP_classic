#include "Monocrystall.h"
#include <cmath>
#define M_PI (3.14159265358979323846264338327950288)
#define M_SQRT2    1.41421356237309504880
#define M_SQRT3    1.73205080756887729352
#define MEGA 1.0e6

Monocrystall::Monocrystall() : Polycrystall(), O({ 3, 3 }) {}


Tensor Monocrystall::to_KSK(const Tensor& T) {
    if (T.dimensions != std::vector<size_t>{3, 3} || O.dimensions != std::vector<size_t>{3, 3}) {
        throw std::invalid_argument("Both T and O must be 3x3 tensors.");
    }

    Tensor T_new({ 3, 3 });
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            double value = 0.0;
            for (size_t k = 0; k < 3; ++k) {
                for (size_t l = 0; l < 3; ++l) {
                    value += O.get({ k, i }) * O.get({ l, j }) * T.get({ k, l });
                }
            }
            T_new.set({ i, j }, value);
        }
    }
    return T_new;
}

Tensor Monocrystall::from_KSK(const Tensor& T)
{
    if (T.dimensions != std::vector<size_t>{3, 3} || O.dimensions != std::vector<size_t>{3, 3}) {
        throw std::invalid_argument("Both T and O must be 3x3 tensors.");
    }

    Tensor T_new({ 3, 3 });
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            double value = 0.0;
            for (size_t k = 0; k < 3; ++k) {
                for (size_t l = 0; l < 3; ++l) {
                    value += O.get({ i, k }) * O.get({ j, l }) * T.get({ k, l });
                }
            }
            T_new.set({ i, j }, value);
        }
    }
    return T_new;
}

void Monocrystall::uniform_distribution() {
    double a = random_double(0.0, 2.0 * M_PI);
    double b = acos(random_double(-1.0, 1.0));
    double g = random_double(0.0, 2.0 * M_PI);
    double ca = cos(a);
    double sa = sin(a);
    double cb = cos(b);
    double sb = sin(b);
    double cg = cos(g);
    double sg = sin(g);

    O.set({ 0, 0 }, ca * cb * cg - sa * sg);
    O.set({ 0, 1 }, -cg * sa - ca * cb * sg);
    O.set({ 0, 2 }, ca * sb);
    O.set({ 1, 0 }, cb * cg * sa + ca * sg);
    O.set({ 1, 1 }, ca * cg - cb * sa * sg);
    O.set({ 1, 2 }, sa * sb);
    O.set({ 2, 0 }, -cg * sb);
    O.set({ 2, 1 }, sb * sg);
    O.set({ 2, 2 }, cb);
}

double Monocrystall::random_double(double min, double max) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(min, max);
    return dis(gen);
}

std::vector<double> Monocrystall::from_KSK(const std::vector<double>& vec) {
    if (vec.size() != 3) {
        throw std::invalid_argument("Input vector must be of size 3.");
    }

    std::vector<double> result(3, 0.0);
    for (size_t i = 0; i < 3; i++) {
        for (size_t j = 0; j < 3; j++) {
            result[i] += O.get({ i, j }) * vec[j];
        }
    }
    return result;
}

std::vector<double> Monocrystall::to_KSK(const std::vector<double>& vec) {
    if (vec.size() != 3) {
        throw std::invalid_argument("Input vector must be of size 3.");
    }

    std::vector<double> result(3, 0.0);
    for (size_t i = 0; i < 3; i++) {
        for (size_t j = 0; j < 3; j++) {
            result[i] += O.get({ j, i }) * vec[j];
        }
    }
    return result;
}

void Monocrystall::fillElast4D(const Params& param) {
    double koef = param.get_param("koef");
    double c11 = koef * param.get_param("c11") / MEGA;
    double c12 = koef * param.get_param("c12") / MEGA;
    double c44 = koef * param.get_param("c44") / MEGA;

    for (size_t i = 0; i < 3; i++) {
        for (size_t j = 0; j < 3; j++) {
            for (size_t k = 0; k < 3; k++) {
                for (size_t s = 0; s < 3; s++) {
                    if ((i == j) && (j == k) && (k == s)) {
                        elast4D.set({ i, j, k, s }, c11);
                    }
                    else if ((i == j) && (k == s)) {
                        elast4D.set({ i, j, k, s }, c12);
                    }
                    else if (((i == k) && (j == s)) || ((i == s) && (j == k))) {
                        elast4D.set({ i, j, k, s }, c44);
                    }
                    else {
                        elast4D.set({ i, j, k, s }, 0.0);
                    }
                }
            }
        }
    }

    // Convert elast4D to elast2D
    elast_4D_to_2D();
}