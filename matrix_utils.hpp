#ifndef MATRIX_UTILS_HPP
#define MATRIX_UTILS_HPP

#include <cstddef>
#include <stdexcept>

// Estructura para manejar matrices
struct Matrix {
    double** data;
    size_t rows;
    size_t cols;
};

// Crear matriz con ceros
Matrix create_matrix(size_t rows, size_t cols) {
    Matrix m;
    m.rows = rows;
    m.cols = cols;
    m.data = new double*[rows];
    for (size_t i = 0; i < rows; i++) {
        m.data[i] = new double[cols]();
    }
    return m;
}

// Liberar memoria
void free_matrix(Matrix& m) {
    for (size_t i = 0; i < m.rows; i++) {
        delete[] m.data[i];
    }
    delete[] m.data;
    m.data = nullptr;
}

// Copiar matriz
Matrix copy_matrix(const Matrix& src) {
    Matrix dst = create_matrix(src.rows, src.cols);
    for (size_t i = 0; i < src.rows; i++) {
        for (size_t j = 0; j < src.cols; j++) {
            dst.data[i][j] = src.data[i][j];
        }
    }
    return dst;
}

// Sumar dos matrices
Matrix add_matrix(const Matrix& A, const Matrix& B) {
    if (A.rows != B.rows || A.cols != B.cols) {
        throw std::invalid_argument("Dimensiones incompatibles para suma");
    }
    Matrix C = create_matrix(A.rows, A.cols);
    for (size_t i = 0; i < A.rows; i++) {
        for (size_t j = 0; j < A.cols; j++) {
            C.data[i][j] = A.data[i][j] + B.data[i][j];
        }
    }
    return C;
}

// Restar dos matrices
Matrix subtract_matrix(const Matrix& A, const Matrix& B) {
    if (A.rows != B.rows || A.cols != B.cols) {
        throw std::invalid_argument("Dimensiones incompatibles para resta");
    }
    Matrix C = create_matrix(A.rows, A.cols);
    for (size_t i = 0; i < A.rows; i++) {
        for (size_t j = 0; j < A.cols; j++) {
            C.data[i][j] = A.data[i][j] - B.data[i][j];
        }
    }
    return C;
}

#endif