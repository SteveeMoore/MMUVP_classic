#include <iostream>
#include <vector>
#include "Tensor.h"
#include "Params.h"
#include "Polycrystall.h"
#include "Monocrystall.h"

int main() {
    Params params("input\\params.input");

    Polycrystall polycrystall;

    // Initialize the Monocrystall object
    std::vector<Monocrystall> monocrystalls(1000);

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
    

    return 0;
}
