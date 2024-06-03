#include <iostream>
#include "Tensor.h"
#include "Params.h"
#include "Polycrystall.h"
#include "Monocrystall.h"

int main() {
    Params params("input\\params.input");

    Polycrystall polycrystall;

    // Initialize the Monocrystall object
    Monocrystall monocrystall;

    // Fill elast4D and convert to elast2D
    monocrystall.fillElast4D(params);
    
    // Print the resulting elast2D tensor
    std::cout << "Converted elast2D Tensor:" << std::endl;
    monocrystall.elast2D.print();

    return 0;
}
