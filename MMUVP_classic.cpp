#include <iostream>
#include <vector>
#include "Tensor.h"
#include "Params.h"
#include "Polycrystall.h"
#include "Monocrystall.h"
#include "Common.h"

int main() {
    Common::clear_output_folder();
    Params params("input\\params.input");

    Polycrystall polycrystall;

    // Initialize the Monocrystall object
    std::vector<Monocrystall> monocrystalls((int)params.get_param("grain_num"));

    //Инициализация
    for (auto& mono: monocrystalls)
    {
        mono.fillElast4D(params);
        mono.elast2D.save_to_file("output\\P.dat");
        mono.uniform_distribution();
        mono.O.save_to_file("output\\O.dat");
    }
    polycrystall.calcAverageElast4D(monocrystalls);
    polycrystall.elast2D.save_to_file("output\\LSK_P.dat");
    polycrystall.true_uniaxial(params.get_param("L11"));

    for (auto& mono : monocrystalls)
    {
        mono.D = mono.to_KSK(polycrystall.D);
        mono.D.save_to_file("output\\D_mono.dat");
    }

    polycrystall.save_SE11_to_file("output\\SE11.dat", polycrystall.S_mean);
    polycrystall.save_SE11_to_file("output\\SE11_poly.dat", polycrystall.S);

    double dt = params.get_param("dt");
    size_t step = 0;
    size_t write_step = (size_t)params.get_param("write_step");
    while (polycrystall.E.get({ 0,0 }) < 0.005)
    {
        step += 1;
        for (auto& mono : monocrystalls)
        {
            mono.calcHookeLaw4D();
            mono.integrateS_dot(dt);
            mono.integrateD(dt);
        }

        polycrystall.calcHookeLaw4D();
        polycrystall.calcAverageS(monocrystalls);
        polycrystall.integrateS_dot(dt);
        polycrystall.integrateD(dt);
        if (step % 3 == 0) {
            polycrystall.save_SE11_to_file("output\\SE11.dat", polycrystall.S_mean);
            polycrystall.save_SE11_to_file("output\\SE11_poly.dat", polycrystall.S);
            polycrystall.S.save_to_file("output\\S_poly.dat");
            polycrystall.E.print();
        }
    }

    return 0;
}
