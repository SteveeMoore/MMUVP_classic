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

void Polycrystall::calcHookeLaw() {
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

