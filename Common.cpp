#include "Common.h"

void Common::clear_output_folder(const std::string& directoryPath) {
    //try {
        // ������� ����� � �� ����������
    std::filesystem::remove_all(directoryPath);
    //std::cout << "����� ������� �������: " << directoryPath << std::endl;

    // ������� ����� ������
    std::filesystem::create_directory(directoryPath);
    //std::cout << "����� ������� ������� ������: " << directoryPath << std::endl;
//}
//catch (const std::exception& e) {
//    std::cerr << "������: " << e.what() << std::endl;
//}
}