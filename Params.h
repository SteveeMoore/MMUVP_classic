#pragma once

#include <iostream>
#include <fstream>
#include <unordered_map>
#include <string>

class Params {
public:
    std::unordered_map<std::string, double> value;

    Params(const std::string& filename);

    // Method to get a parameter value by key
    double get_param(const std::string& key) const;

    // Method to set a parameter value by key
    void set_param(const std::string& key, double val);

private:
    std::unordered_map<std::string, double> from_file_to_map(const std::string& filename);
};

