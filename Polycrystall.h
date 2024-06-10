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

    /*---------�������������---------*/
    //������� ������� ������� ������� � ����� 6x6
    void elast_4D_to_2D();
    

    
    /*---------����������---------*/
    //���������� ������ ���� � ����������� �����
    void calcHookeLaw2D();
    //���������� ������ ���� � ������ �����
    void calcHookeLaw4D();

    /*---------��������������---------*/
    //�������������� ������� �������� ����������
    void integrateS_dot(double dt);
    //�������������� ������� �������� ����������
    void integrateD(double dt);

    /*---------����������---------*/
    //���������� ������� ������� �������
    void calcAverageElast4D(std::vector<Monocrystall>& monocrystals);
    //���������� ���������� �� ������������. 
    void calcAverageS(std::vector<Monocrystall>& monocrystals);

    /*---------������� ������� ��������� � ��������� �������---------*/
    void true_uniaxial(double L11);


    /*---------����� � ����---------*/
    void save_SE11_to_file(const std::string& filename, Tensor& tensor);
    void save_intensity_to_file(const std::string& filename, Tensor& tensor);


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
};
