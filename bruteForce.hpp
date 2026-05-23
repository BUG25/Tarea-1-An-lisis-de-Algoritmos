#ifndef BRUTEFORCE_HPP
#define BRUTEFORCE_HPP

#include "matrix_utils.hpp"

// Multiplicación de matrices por fuerza bruta (algoritmo cúbico)
Matrix brute_force(const Matrix& A, const Matrix& B) {
    if (A.cols != B.rows) {
        throw std::invalid_argument("Dimensiones incompatibles para multiplicación");
    }
    
    size_t n1 = A.rows;    // filas de A
    size_t n2 = A.cols;    // columnas de A = filas de B
    size_t n3 = B.cols;    // columnas de B
    
    Matrix C = create_matrix(n1, n3);
    
    for (size_t i = 0; i < n1; i++) {
        for (size_t j = 0; j < n2; j++) {
            for (size_t k = 0; k < n3; k++) {
                C.data[i][k] += A.data[i][j] * B.data[j][k];
            }
        }
    }
    
    return C;
}

#endif