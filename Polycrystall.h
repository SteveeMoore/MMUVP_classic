#pragma once
#include <string>
#include "Tensor.h"
#include "Params.h"
//#include "Monocrystall.h"

class Monocrystall;

class Polycrystall {
public:
    Polycrystall();

    /*---------�������������---------*/
    //������� ������� ������� ������� � ����� 6x6
    void elast_4D_to_2D();
    //���������� ������� ������� �������
    void calcAverageElast4D(std::vector<Monocrystall>& monocrystals);
    

    
    /*---------����������---------*/
    //���������� ������ ���� � ����������� �����
    void calcHookeLaw2D();
    //���������� ������ ���� � ������ �����
    void calcHookeLaw4D();

    /*---------��������������---------*/
    //�������������� ������� �������� ����������
    void integrateS_dot(double dt);
    //�������������� ������� �������� ����������
    void integrateL(double dt);


    /*---------������� ������� ��������� � ��������� �������---------*/
    void true_uniaxial(double L11);



public:
    //������ ���������� ��������
    Tensor L;
    //������������ ����� ���������� ��������
    SymTensor D;
    //����������� ����������
    SymTensor E;

    //������ ����������
    SymTensor S;
    //������ �������� ��������� ����������
    SymTensor S_dot;
    //����������� ������ ����������
    SymTensor S_mean;
    //������ ������� ������� 3�3�3�3
    Tensor elast4D;
    //������ ������� ������� 6�6
    Tensor elast2D;

protected:
    std::vector<double> solveGaussStrictOrder(const std::vector<std::vector<double>>& matrix, const std::vector<double>& results);
};
