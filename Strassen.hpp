#ifndef STRASSEN_HPP
#define STRASSEN_HPP

#include "BruteForce.hpp"

// Dividir matriz en 4 submatrices
void split_matrix(const Matrix& M, Matrix& a, Matrix& b, Matrix& c, Matrix& d) {
    size_t half = M.rows / 2;
    
    a = create_matrix(half, half);
    b = create_matrix(half, half);
    c = create_matrix(half, half);
    d = create_matrix(half, half);
    
    for (size_t i = 0; i < half; i++) {
        for (size_t j = 0; j < half; j++) {
            a.data[i][j] = M.data[i][j];
            b.data[i][j] = M.data[i][j + half];
            c.data[i][j] = M.data[i + half][j];
            d.data[i][j] = M.data[i + half][j + half];
        }
    }
}

// Unir 4 submatrices en una matriz
Matrix join_matrices(const Matrix& a, const Matrix& b, const Matrix& c, const Matrix& d) {
    size_t half = a.rows;
    size_t full = half * 2;
    Matrix M = create_matrix(full, full);
    
    for (size_t i = 0; i < half; i++) {
        for (size_t j = 0; j < half; j++) {
            M.data[i][j] = a.data[i][j];
            M.data[i][j + half] = b.data[i][j];
            M.data[i + half][j] = c.data[i][j];
            M.data[i + half][j + half] = d.data[i][j];
        }
    }
    
    return M;
}

// Algoritmo de Strassen
Matrix strassen(const Matrix& A, const Matrix& B) {
    if (A.rows == 1) {
        return brute_force(A, B);
    }
    
    // Dividir matrices A y B
    Matrix a, b, c, d;
    Matrix e, f, g, h;
    split_matrix(A, a, b, c, d);
    split_matrix(B, e, f, g, h);
    
    // Calcular los 7 productos de Strassen
    // p1 = (a+d)(e+h)
    Matrix a_plus_d = add_matrix(a, d);
    Matrix e_plus_h = add_matrix(e, h);
    Matrix p1 = strassen(a_plus_d, e_plus_h);
    free_matrix(a_plus_d);
    free_matrix(e_plus_h);
    
    // p2 = d(g-e)
    Matrix g_minus_e = subtract_matrix(g, e);
    Matrix p2 = strassen(d, g_minus_e);
    free_matrix(g_minus_e);
    
    // p3 = (a+b)h
    Matrix a_plus_b = add_matrix(a, b);
    Matrix p3 = strassen(a_plus_b, h);
    free_matrix(a_plus_b);
    
    // p4 = (b-d)(g+h)
    Matrix b_minus_d = subtract_matrix(b, d);
    Matrix g_plus_h = add_matrix(g, h);
    Matrix p4 = strassen(b_minus_d, g_plus_h);
    free_matrix(b_minus_d);
    free_matrix(g_plus_h);
    
    // p5 = a(f-h)
    Matrix f_minus_h = subtract_matrix(f, h);
    Matrix p5 = strassen(a, f_minus_h);
    free_matrix(f_minus_h);
    
    // p6 = (c+d)e
    Matrix c_plus_d = add_matrix(c, d);
    Matrix p6 = strassen(c_plus_d, e);
    free_matrix(c_plus_d);
    
    // p7 = (a-c)(e+f)
    Matrix a_minus_c = subtract_matrix(a, c);
    Matrix e_plus_f = add_matrix(e, f);
    Matrix p7 = strassen(a_minus_c, e_plus_f);
    free_matrix(a_minus_c);
    free_matrix(e_plus_f);
    
    // Calcular las 4 submatrices resultado
    // C11 = p1 + p2 - p3 + p4
    Matrix temp1 = add_matrix(p1, p2);
    Matrix temp2 = subtract_matrix(temp1, p3);
    Matrix C11 = add_matrix(temp2, p4);
    free_matrix(temp1);
    free_matrix(temp2);
    
    // C12 = p5 + p3
    Matrix C12 = add_matrix(p5, p3);
    
    // C21 = p6 + p2
    Matrix C21 = add_matrix(p6, p2);
    
    // C22 = p5 + p1 - p6 - p7
    Matrix temp3 = add_matrix(p5, p1);
    Matrix temp4 = subtract_matrix(temp3, p6);
    Matrix C22 = subtract_matrix(temp4, p7);
    free_matrix(temp3);
    free_matrix(temp4);
    
    // Unir las 4 submatrices
    Matrix C = join_matrices(C11, C12, C21, C22);
    
    // Liberar memoria de submatrices y productos intermedios
    free_matrix(a); free_matrix(b); free_matrix(c); free_matrix(d);
    free_matrix(e); free_matrix(f); free_matrix(g); free_matrix(h);
    free_matrix(p1); free_matrix(p2); free_matrix(p3); free_matrix(p4);
    free_matrix(p5); free_matrix(p6); free_matrix(p7);
    free_matrix(C11); free_matrix(C12); free_matrix(C21); free_matrix(C22);
    
    return C;
}

#endif