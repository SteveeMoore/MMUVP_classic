#include "Tensor.h"

Tensor::Tensor() : total_size(0) {}

Tensor::Tensor(const std::vector<size_t>& dims) : dimensions(dims) {
    total_size = 1;
    for (size_t dim : dims) {
        total_size *= dim;
    }
    data.resize(total_size, 0);
}

size_t Tensor::get_index(const std::vector<size_t>& indices) const {
    size_t index = 0;
    size_t multiplier = 1;
    for (int i = dimensions.size() - 1; i >= 0; --i) {
        index += indices[i] * multiplier;
        multiplier *= dimensions[i];
    }
    return index;
}

double Tensor::get(const std::vector<size_t>& indices) const {
    if (indices.size() != dimensions.size()) {
        throw std::invalid_argument("Incorrect number of indices.");
    }
    return data[get_index(indices)];
}

void Tensor::set(const std::vector<size_t>& indices, double value) {
    if (indices.size() != dimensions.size()) {
        throw std::invalid_argument("Incorrect number of indices.");
    }
    data[get_index(indices)] = value;
}

Tensor Tensor::operator+(const Tensor& other) const {
    if (dimensions != other.dimensions) {
        throw std::invalid_argument("Tensors must have the same dimensions.");
    }
    Tensor result(dimensions);
    for (size_t i = 0; i < total_size; ++i) {
        result.data[i] = data[i] + other.data[i];
    }
    return result;
}

Tensor Tensor::operator*(double scalar) const {
    Tensor result(dimensions);
    for (size_t i = 0; i < total_size; ++i) {
        result.data[i] = data[i] * scalar;
    }
    return result;
}

Tensor Tensor::operator+(double scalar) const {
    Tensor result(dimensions);
    for (size_t i = 0; i < total_size; ++i) {
        result.data[i] = data[i] + scalar;
    }
    return result;
}

void Tensor::print() const {
    print_recursive(0, {});
}

void Tensor::print_recursive(size_t dim, std::vector<size_t> indices) const {
    if (dim == dimensions.size()) {
        std::cout << get(indices) << " ";
        return;
    }

    std::cout << "[ ";
    for (size_t i = 0; i < dimensions[dim]; ++i) {
        indices.push_back(i);
        print_recursive(dim + 1, indices);
        indices.pop_back();
    }
    std::cout << "] ";
    if (dim == 0) {
        std::cout << std::endl;
    }
}

void Tensor::save_to_file(const std::string& filename) const {
    std::ofstream file;
    file.open(filename, std::ios_base::app); // Open file in append mode
    if (!file) {
        throw std::runtime_error("Could not open file.");
    }
    save_to_file_recursive(file, 0, {});
    file << "\n"; // Add a newline at the end of the tensor
    file.close();
}

void Tensor::save_to_file_recursive(std::ofstream& file, size_t dim, std::vector<size_t> indices) const {
    if (dim == dimensions.size()) {
        file << get(indices) << " ";
        return;
    }

    for (size_t i = 0; i < dimensions[dim]; ++i) {
        indices.push_back(i);
        save_to_file_recursive(file, dim + 1, indices);
        indices.pop_back();
    }
}

// Symmetric Tensor Implementation

SymTensor::SymTensor(const std::vector<size_t>& dims) : Tensor(dims) {
    if (dims.size() != 2 || dims[0] != dims[1]) {
        throw std::invalid_argument("Symmetric tensor must be a square matrix.");
    }
}

SymTensor::SymTensor(const Tensor& other) : Tensor(other) {
    if (dimensions.size() != 2 || dimensions[0] != dimensions[1]) {
        throw std::invalid_argument("Symmetric tensor must be a square matrix.");
    }
    symmetrize();
}

void SymTensor::symmetrize() {
    for (size_t i = 0; i < dimensions[0]; ++i) {
        for (size_t j = i + 1; j < dimensions[1]; ++j) {
            double avg = (get({ i, j }) + get({ j, i })) / 2.0;
            set({ i, j }, avg);
            set({ j, i }, avg);
        }
    }
}

std::vector<double> SymTensor::to_vector() const {
    std::vector<double> result;
    result.push_back(get({ 0, 0 }));
    result.push_back(get({ 1, 1 }));
    result.push_back(get({ 2, 2 }));
    result.push_back(get({ 0, 1 }));
    result.push_back(get({ 0, 2 }));
    result.push_back(get({ 1, 2 }));
    return result;
}
