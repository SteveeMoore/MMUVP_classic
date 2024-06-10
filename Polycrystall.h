#pragma once
#include <string>
#include "Tensor.h"
#include "Params.h"
#include "Solver.h"
//#include "Monocrystall.h"

class Monocrystall;

class Polycrystall {
public:
    Polycrystall();

    /*---------Инициализация---------*/
    //Перевод тензора упругих модулей в форму 6x6
    void elast_4D_to_2D();
    

    
    /*---------Вычисления---------*/
    //Вычисление закона Гука в сокращенной форме
    void calcHookeLaw2D();
    //Вычисление закона Гука в полной форме
    void calcHookeLaw4D();

    /*---------Интегрирование---------*/
    //Интегрирование тензора скорости напряжений
    void integrateS_dot(double dt);
    //Интегрирование тензора скорости деформаций
    void integrateD(double dt);

    /*---------Осреднение---------*/
    //Осреднение тензора упругих свойств
    void calcAverageElast4D(std::vector<Monocrystall>& monocrystals);
    //Осреднение напряжений по кристаллитам. 
    void calcAverageS(std::vector<Monocrystall>& monocrystals);

    /*---------Функции задания начальных и граничных условий---------*/
    void true_uniaxial(double L11);


    /*---------Вывод в файл---------*/
    void save_SE11_to_file(const std::string& filename, Tensor& tensor);
    void save_intensity_to_file(const std::string& filename, Tensor& tensor);


public:
    //Тензор деформации скорости
    Tensor L;
    //Симметричная часть деформации скорости
    SymTensor D;
    //Накопленные деформации
    SymTensor E;

    //Тензор напряжений
    SymTensor S;
    //Тензор скорости изменения напряжений
    SymTensor S_dot;
    //Осредненный тензор напряжений
    SymTensor S_mean;
    //Тензор упругих модулей 3х3х3х3
    Tensor elast4D;
    //Тензор упругих модулей 6х6
    Tensor elast2D;
};
