#ifndef HYBRID_HPP
#define HYBRID_HPP

#include "bruteforce.hpp"
#include "strassen.hpp"

Matrix hybrido(const Matrix& A, const Matrix& B, size_t n0) {
    
    if (A.rows <= n0) {
        return brute_force(A, B);
    } else {
        return strassen(A, B, n0);
    }
}


#endif 