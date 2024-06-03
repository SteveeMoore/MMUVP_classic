#pragma once

#include <iostream>
#include <vector>
#include <stdexcept>
#include <fstream>

class Tensor {
protected:
    std::vector<double> data;
    
    size_t total_size;

    size_t get_index(const std::vector<size_t>& indices) const;

public:
    Tensor();
    Tensor(const std::vector<size_t>& dims);

    double get(const std::vector<size_t>& indices) const;
    void set(const std::vector<size_t>& indices, double value);

    Tensor operator+(const Tensor& other) const;
    Tensor operator*(double scalar) const;
    Tensor operator+(double scalar) const;

    void print() const;
    void save_to_file(const std::string& filename) const;

    std::vector<size_t> dimensions;
protected:
    void print_recursive(size_t dim, std::vector<size_t> indices) const;
    void save_to_file_recursive(std::ofstream& file, size_t dim, std::vector<size_t> indices) const;
};


class SymTensor : public Tensor {
public:
    SymTensor(const std::vector<size_t>& dims);
    SymTensor(const Tensor& other);

    //TODO: Сделать универсальнее
    void symmetrize();
    std::vector<double> to_vector() const;
};

