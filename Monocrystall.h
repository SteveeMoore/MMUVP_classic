#pragma once
#include "Polycrystall.h"

#include <random>

class Monocrystall : public Polycrystall {
public:
    Monocrystall();

    //�������������
    void uniform_distribution();


    //������� � ���
    Tensor to_KSK(const Tensor& T);
    std::vector<double> to_KSK(const std::vector<double>& vec);
    //������� � ���
    Tensor from_KSK(const Tensor& T);
    std::vector<double> from_KSK(const std::vector<double>& vec);

    void fillElast4D(const Params& param);
    
    
private:
    static double random_double(double min, double max);
public:
    // ������ ���������� �
    Tensor O;
};
