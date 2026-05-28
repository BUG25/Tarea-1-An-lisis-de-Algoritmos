#include <iostream>
#include "bruteforce.hpp"
#include "Strassen.hpp"

void print_matrix(const Matrix& M) {
    for (size_t i = 0; i < M.rows; i++) {
        for (size_t j = 0; j < M.cols; j++) {
            std::cout << M.data[i][j] << "\t";
        }
        std::cout << std::endl;
    }
}

void fill_random(Matrix& M) {
    for (size_t i = 0; i < M.rows; i++) {
        for (size_t j = 0; j < M.cols; j++) {
            M.data[i][j] = rand() % 10;
        }
    }
}

int main() {
    int e = 2;
    for (int t = 1; t <= 3; t++) {

         Matrix A = create_matrix(e,e);
        Matrix B = create_matrix(e,e);

        fill_random(A);
        fill_random(B);

        Matrix C = brute_force(A, B);
        Matrix D = strassen(A, B, 1);

        std::cout << "\n====================\n";
        std::cout << "Instancia " << t << "\n";
        std::cout << "====================\n";

        std::cout << "\nMatriz A:\n";
        print_matrix(A);

        std::cout << "\nMatriz B:\n";
        print_matrix(B);

        std::cout << "\nBrute Force:\n";
        print_matrix(C);

        std::cout << "\nStrassen:\n";
        print_matrix(D);

        free_matrix(A);
        free_matrix(B);
        free_matrix(C);
        free_matrix(D);
        e=e*2;
    }

    return 0;
}