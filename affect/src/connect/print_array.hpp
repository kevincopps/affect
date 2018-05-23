#include <iostream>

template<class T>
inline void print_1d_array(const char* name, const T* array, size_t n) {
    std::cout << std::endl << name << "(" << n << "):" << std::endl << "    [";
    for (size_t i = 0; i < n; ++i) {
        std::cout << array[i] << ", ";
    }
    std::cout << "]" << std::endl;
}

template<>
inline void print_1d_array<int8_t>(const char* name, const int8_t* array, size_t n) {
    std::cout << std::endl << name << "(" << n << "):" << std::endl << "    [";
    for (size_t i = 0; i < n; ++i) {
        std::cout << (int)array[i] << ", ";
    }
    std::cout << "]" << std::endl;
}


template<class T>
inline void print_2d_array(const char* name, const T* array, size_t m, size_t n) {
    std::cout << std::endl << name << "(" << m << ", " << n << "):" << std::endl << "   [" << std::endl;
    for (size_t i = 0; i < m; ++i) {
        std::cout << "    [";
        size_t ioffset = i * m;
        for (size_t j = 0; j < n; ++j) {
            std::cout << array[ioffset + j] << ", ";
        }
        std::cout << "]," << std::endl;
    }
    std::cout << "   ]" << std::endl;
}