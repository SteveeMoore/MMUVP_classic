#include "Common.h"

void Common::clear_output_folder(const std::string& directoryPath) {
    //try {
        // Удаляем папку и ее содержимое
    std::filesystem::remove_all(directoryPath);
    //std::cout << "Папка успешно удалена: " << directoryPath << std::endl;

    // Создаем папку заново
    std::filesystem::create_directory(directoryPath);
    //std::cout << "Папка успешно создана заново: " << directoryPath << std::endl;
//}
//catch (const std::exception& e) {
//    std::cerr << "Ошибка: " << e.what() << std::endl;
//}
}