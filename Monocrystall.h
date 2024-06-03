#pragma once
#include "Polycrystall.h"

#include <random>

class Monocrystall : public Polycrystall {
public:
    Monocrystall();

    //Инициализация
    void uniform_distribution();


    //Перевод в КСК
    Tensor to_KSK(const Tensor& T);
    std::vector<double> to_KSK(const std::vector<double>& vec);
    //Перевод в ЛСК
    Tensor from_KSK(const Tensor& T);
    std::vector<double> from_KSK(const std::vector<double>& vec);

    void fillElast4D(const Params& param);
    
    
private:
    static double random_double(double min, double max);
public:
    // Тензор ориентации О
    Tensor O;
};
