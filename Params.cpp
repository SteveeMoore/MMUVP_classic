#include "Params.h"

Params::Params(const std::string& filename) {
    value = from_file_to_map("input/params.input");
}

std::unordered_map<std::string, double> Params::from_file_to_map(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Could not open the file: " << filename << std::endl;
        return {};
    }

    std::unordered_map<std::string, double> data;
    std::string line;

    while (std::getline(file, line)) {
        size_t quote_pos1 = line.find('"');
        size_t quote_pos2 = line.find('"', quote_pos1 + 1);
        if (quote_pos1 == std::string::npos || quote_pos2 == std::string::npos)
            continue; // Skip lines without quotes

        std::string key = line.substr(quote_pos1 + 1, quote_pos2 - quote_pos1 - 1); // Extracting key
        double value = std::stod(line.substr(quote_pos2 + 3)); // Extracting value after the colon
        data[key] = value;
    }
    file.close();
    return data;
}

double Params::get_param(const std::string& key) const {
    auto it = value.find(key);
    if (it != value.end()) {
        return it->second;
    }
    throw std::invalid_argument("Key not found: " + key);
}

void Params::set_param(const std::string& key, double val) {
    value[key] = val;
}
