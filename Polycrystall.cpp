#include "Polycrystall.h"
#include "Monocrystall.h"

Polycrystall::Polycrystall()
    : L({ 3, 3 }),
    D({ 3, 3 }),
    E({ 3, 3 }),
    S({ 3, 3 }),
    S_dot({ 3, 3 }),
    S_mean({ 3, 3 }),
    elast4D({ 3, 3, 3, 3 }),
    elast2D({ 6, 6 }) 
{}

void Polycrystall::elast_4D_to_2D() {
    std::vector<std::vector<size_t>> indices = {
        {0, 0, 0, 0}, {0, 0, 1, 1}, {0, 0, 2, 2}, {0, 0, 1, 2}, {0, 0, 2, 0}, {0, 0, 0, 1},
        {1, 1, 0, 0}, {1, 1, 1, 1}, {1, 1, 2, 2}, {1, 1, 1, 2}, {1, 1, 2, 0}, {1, 1, 0, 1},
        {2, 2, 0, 0}, {2, 2, 1, 1}, {2, 2, 2, 2}, {2, 2, 1, 2}, {2, 2, 2, 0}, {2, 2, 0, 1},
        {1, 2, 0, 0}, {1, 2, 1, 1}, {1, 2, 2, 2}, {1, 2, 1, 2}, {1, 2, 2, 0}, {1, 2, 0, 1},
        {2, 0, 0, 0}, {2, 0, 1, 1}, {2, 0, 2, 2}, {2, 0, 1, 2}, {2, 0, 2, 0}, {2, 0, 0, 1},
        {0, 1, 0, 0}, {0, 1, 1, 1}, {0, 1, 2, 2}, {0, 1, 1, 2}, {0, 1, 2, 0}, {0, 1, 0, 1}
    };

    for (size_t i = 0; i < 6; ++i) {
        for (size_t j = 0; j < 6; ++j) {
            size_t idx = i * 6 + j;
            elast2D.set({ i, j }, elast4D.get({ indices[idx][0], indices[idx][1], indices[idx][2], indices[idx][3] }));
        }
    }
}

void Polycrystall::calcAverageElast4D(std::vector<Monocrystall>& monocrystals) {
    Tensor averageElast4D({ 3, 3, 3, 3 });
    size_t numCrystals = monocrystals.size();

    if (numCrystals == 0) {
        throw std::invalid_argument("The vector of monocrystals is empty.");
    }

    for (auto& mono : monocrystals) {
        auto mono_elast = (mono.elast_from_KSK());
        for (size_t i = 0; i < 3; i++) {
            for (size_t j = 0; j < 3; j++) {
                for (size_t k = 0; k < 3; k++) {
                    for (size_t s = 0; s < 3; s++) {
                        double currentValue = averageElast4D.get({ i, j, k, s });
                        double monoValue = mono_elast.get({ i, j, k, s });
                        averageElast4D.set({ i, j, k, s }, currentValue + monoValue);
                    }
                }
            }
        }
    }

    for (size_t i = 0; i < 3; i++) {
        for (size_t j = 0; j < 3; j++) {
            for (size_t k = 0; k < 3; k++) {
                for (size_t s = 0; s < 3; s++) {
                    double currentValue = averageElast4D.get({ i, j, k, s });
                    elast4D.set({ i, j, k, s }, currentValue / (double)numCrystals);
                }
            }
        }
    }

    elast_4D_to_2D();
}

void Polycrystall::calcHookeLaw2D() {
    std::vector<double> sigma_rate_vector(6, 0.0);
    std::vector<double> strain_rate_vector = D.to_vector();

    for (size_t i = 0; i < 6; i++) {
        for (size_t j = 0; j < 6; j++) {
            sigma_rate_vector[i] += elast2D.get({ i, j }) * strain_rate_vector[j];
        }
    }

    S_dot.set({ 0, 0 }, sigma_rate_vector[0]);
    S_dot.set({ 1, 1 }, sigma_rate_vector[1]);
    S_dot.set({ 2, 2 }, sigma_rate_vector[2]);
    S_dot.set({ 0, 1 }, sigma_rate_vector[3]);
    S_dot.set({ 1, 0 }, sigma_rate_vector[3]);
    S_dot.set({ 1, 2 }, sigma_rate_vector[5]);
    S_dot.set({ 2, 1 }, sigma_rate_vector[5]);
    S_dot.set({ 2, 0 }, sigma_rate_vector[4]);
    S_dot.set({ 0, 2 }, sigma_rate_vector[4]);
}

void Polycrystall::calcHookeLaw4D()
{
    // Initialize sigma_rate tensor
    SymTensor sigma_rate({ 3, 3 });
    sigma_rate.set({ 0, 0 }, 0.0);
    sigma_rate.set({ 1, 1 }, 0.0);
    sigma_rate.set({ 2, 2 }, 0.0);
    sigma_rate.set({ 0, 1 }, 0.0);
    sigma_rate.set({ 1, 0 }, 0.0);
    sigma_rate.set({ 1, 2 }, 0.0);
    sigma_rate.set({ 2, 1 }, 0.0);
    sigma_rate.set({ 2, 0 }, 0.0);
    sigma_rate.set({ 0, 2 }, 0.0);

    // Calculate sigma_rate tensor using the elast4D tensor
    for (size_t i = 0; i < 3; i++) {
        for (size_t j = 0; j < 3; j++) {
            for (size_t k = 0; k < 3; k++) {
                for (size_t l = 0; l < 3; l++) {
                    sigma_rate.set({ i, j }, sigma_rate.get({ i, j }) + elast4D.get({ i, j, k, l }) * D.get({ k, l }));
                }
            }
        }
    }

    // Set the values of S_dot tensor
    S_dot.set({ 0, 0 }, sigma_rate.get({ 0, 0 }));
    S_dot.set({ 1, 1 }, sigma_rate.get({ 1, 1 }));
    S_dot.set({ 2, 2 }, sigma_rate.get({ 2, 2 }));
    S_dot.set({ 0, 1 }, sigma_rate.get({ 0, 1 }));
    S_dot.set({ 1, 0 }, sigma_rate.get({ 0, 1 }));
    S_dot.set({ 1, 2 }, sigma_rate.get({ 1, 2 }));
    S_dot.set({ 2, 1 }, sigma_rate.get({ 1, 2 }));
    S_dot.set({ 2, 0 }, sigma_rate.get({ 2, 0 }));
    S_dot.set({ 0, 2 }, sigma_rate.get({ 2, 0 }));
}

void Polycrystall::integrateS_dot(double dt) {
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            S.set({ i, j }, S.get({ i, j }) + S_dot.get({ i, j }) * dt);
        }
    }
}

void Polycrystall::integrateL(double dt)
{
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            E.set({ i, j }, E.get({ i, j }) + L.get({ i, j }) * dt);
        }
    }
}

void Polycrystall::true_uniaxial(double L11)
{
    std::vector<std::vector<size_t>> indices = {
        {0, 0, 0, 0}, {0, 0, 0, 1}, {0, 0, 0, 2}, {0, 0, 1, 0}, {0, 0, 1, 1}, {0, 0, 1, 2}, {0, 0, 2, 0}, {0, 0, 2, 1}, {0, 0, 2, 2},
        {0, 1, 0, 0}, {0, 1, 0, 1}, {0, 1, 0, 2}, {0, 1, 1, 0}, {0, 1, 1, 1}, {0, 1, 1, 2}, {0, 1, 2, 0}, {0, 1, 2, 1}, {0, 1, 2, 2},
        {0, 2, 0, 0}, {0, 2, 0, 1}, {0, 2, 0, 2}, {0, 2, 1, 0}, {0, 2, 1, 1}, {0, 2, 1, 2}, {0, 2, 2, 0}, {0, 2, 2, 1}, {0, 2, 2, 2},
        {1, 0, 0, 0}, {1, 0, 0, 1}, {1, 0, 0, 2}, {1, 0, 1, 0}, {1, 0, 1, 1}, {1, 0, 1, 2}, {1, 0, 2, 0}, {1, 0, 2, 1}, {1, 0, 2, 2},
        {1, 1, 0, 0}, {1, 1, 0, 1}, {1, 1, 0, 2}, {1, 1, 1, 0}, {1, 1, 1, 1}, {1, 1, 1, 2}, {1, 1, 2, 0}, {1, 1, 2, 1}, {1, 1, 2, 2},
        {1, 2, 0, 0}, {1, 2, 0, 1}, {1, 2, 0, 2}, {1, 2, 1, 0}, {1, 2, 1, 1}, {1, 2, 1, 2}, {1, 2, 2, 0}, {1, 2, 2, 1}, {1, 2, 2, 2},
        {2, 0, 0, 0}, {2, 0, 0, 1}, {2, 0, 0, 2}, {2, 0, 1, 0}, {2, 0, 1, 1}, {2, 0, 1, 2}, {2, 0, 2, 0}, {2, 0, 2, 1}, {2, 0, 2, 2},
        {2, 1, 0, 0}, {2, 1, 0, 1}, {2, 1, 0, 2}, {2, 1, 1, 0}, {2, 1, 1, 1}, {2, 1, 1, 2}, {2, 1, 2, 0}, {2, 1, 2, 1}, {2, 1, 2, 2},
        {2, 2, 0, 0}, {2, 2, 0, 1}, {2, 2, 0, 2}, {2, 2, 1, 0}, {2, 2, 1, 1}, {2, 2, 1, 2}, {2, 2, 2, 0}, {2, 2, 2, 1}, {2, 2, 2, 2}
    };

    std::vector<std::vector<double>> L_test(3, std::vector<double>(9, 0.0));
    L_test[0][0] = L11;


    std::vector<std::vector<double>> sigma_test(3, std::vector<double>(3, 0.0));

    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            for (size_t k = 0; k < 3; ++k) {
                for (size_t l = 0; l < 3; ++l) {
                    sigma_test[i][j] -= elast4D.get({ i,j,k,l }) * L_test[k][l];
                }
            }
        }
    }

    std::vector<double> right_side(8, 0.0);
    right_side[0] = sigma_test[0][1];
    right_side[1] = sigma_test[0][2];
    right_side[2] = sigma_test[1][0];
    right_side[3] = sigma_test[1][1];
    right_side[4] = sigma_test[1][2];
    right_side[5] = sigma_test[2][0];
    right_side[6] = sigma_test[2][1];
    right_side[7] = sigma_test[2][2];

    std::vector<std::vector<double>> left_side(8, std::vector<double>(8, 0.0));

    for (size_t i = 1; i < 9; ++i) {
        for (size_t j = 1; j < 9; ++j) {
            size_t idx = i * 9 + j;
            left_side[i-1][j-1] = elast4D.get({ indices[idx][0], indices[idx][1], indices[idx][2], indices[idx][3] });
        }
    }

    std::vector<double> result(8, 0.0);
    result = Solver::solveGaussSimple(left_side, right_side);

    D.set({ 0, 0 }, L11);
    D.set({ 0, 1 }, result[0]);
    D.set({ 0, 2 }, result[1]);
    D.set({ 1, 0 }, result[2]);
    D.set({ 1, 1 }, result[3]);
    D.set({ 1, 2 }, result[4]);
    D.set({ 2, 0 }, result[5]);
    D.set({ 2, 1 }, result[6]);
    D.set({ 2, 2 }, result[7]);
    D.symmetrize();
    D.save_to_file("output\\D.dat");
}

