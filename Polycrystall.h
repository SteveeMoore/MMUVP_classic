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
    //Осреднение тензора упругих свойств
    void calcAverageElast4D(std::vector<Monocrystall>& monocrystals);
    

    
    /*---------Вычисления---------*/
    //Вычисление закона Гука в сокращенной форме
    void calcHookeLaw2D();
    //Вычисление закона Гука в полной форме
    void calcHookeLaw4D();

    /*---------Интегрирование---------*/
    //Интегрирование тензора скорости напряжений
    void integrateS_dot(double dt);
    //Интегрирование тензора скорости деформаций
    void integrateL(double dt);


    /*---------Функции задания начальных и граничных условий---------*/
    void true_uniaxial(double L11);



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
