// * RECORDAR: no se pueden usar fokin librerías de multiplicación de matrices! debemos implementarlo nosotro!!!!!!!
// n es múltiplo de 2
// si n=1, se hace la multiplicación directamente

#ifndef STRASSEN_HPP
#define STRASSEN_HPP

#include "BruteForce.hpp"

// ALGORITMO DE STRASEEN

// dividimos las matrices en matrices más enanas
void split_matrix(const Matrix& M, Matrix& a, Matrix& b, Matrix& c, Matrix& d) {
    size_t half = M.rows / 2;
    
    // dividimos en 4 submatrice
    // creamos matrices vacías
    a = create_matrix(half, half);
    b = create_matrix(half, half);
    c = create_matrix(half, half);
    d = create_matrix(half, half);
    
    //copiamos datos a estos 4 cuartos de matrix
    for (size_t i = 0; i < half; i++) {
        for (size_t j = 0; j < half; j++) {
            a.data[i][j] = M.data[i][j];                  // 1ra matriz, filas de 0 a n//2, columnas de 0 a n//2
            b.data[i][j] = M.data[i][j + half];          // 2da matriz, filas de 0 a n//2, columnas de n//2 a n
            c.data[i][j] = M.data[i + half][j];         // 3ra matriz, filas de n//2 a n, columnas de 0 a n//2
            d.data[i][j] = M.data[i + half][j + half]; // 4ta matriz, filas de n//2 a n, columnas de n//2 a n
        }
    }
}

// Unimos submatrices
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

Matrix strassen(const Matrix& A, const Matrix& B) {
    if (A.rows == 1) {  // si n=1 se hace la multiplicación directamente
        return brute_force(A, B);
    }

    Matrix a, b, c, d;
    Matrix e, f, g, h;
    split_matrix(A, a, b, c, d);
    split_matrix(B, e, f, g, h);
    
    // Podríamos calcular con strassen los siguientes productos:
    // ae = strassen(a,e)
    // bg = strassen(b,g)
    // af = strassen(a,f)
    // bh = strassen(b,h)
    // ce = strassen(c,e)
    // dg = strassen(d,g)
    // cf = strassen(c,f)
    // dh = strassen(d,h)
    // y luego C11, C12, C21, C22    
    // y unirlas en una sola matriz C  
    // Tendríamos 8 llamaMdas a strassen de complejidaMd T(n/2), y cinco operaciones de complejidaMd Theta(n^2) 
    // T(n) = 8T(n/2) + Theta(n^2) => T(n) = Theta(n^3) por teorema maestro
    // El problema es que Brute Force tiene la misma complejidaMd temporal, lo que hace inoca la implementación de algoritmo de Strassen
    // Por tanto, podemos reescribir de manera tal que:
    // 1. ae+bg = ae+ah+de+dh+dg-de-ah-bh+bg+bh-dg-dh = (a+d)(e+h)+d(g-e)-(a+b)h+(b-d)(g+h)
    // 2. af+bh = af-ah+ah+bh = a(f-h)+h(a+b)
    // 3. ce+dg = ce+de+dg-de = e(c+d)+d(g-e)
    // 4. cf+dh = af-ah+ae+ah+de+dh-ce-de-ae-af+ce+cf = a(f-h)+(a+d)(e+h)-e(c+d)-(a-c)(e+f)
    // Con estas expresiones, podemos definir:
    
    // p1 = (a+d)(e+h)
    Matrix aMd = add_matrix(a, d);  // a+d
    Matrix eMh = add_matrix(e, h);  // e+h
    Matrix p1 = strassen(aMd, eMh);
    free_matrix(aMd);
    free_matrix(eMh);
    
    // p2 = d(g-e)
    Matrix gme = subtract_matrix(g, e);  // g-e
    Matrix p2 = strassen(d, gme);
    free_matrix(gme);
    
    // p3 = (a+b)h
    Matrix aMb = add_matrix(a, b);  // a+b
    Matrix p3 = strassen(aMb, h);
    free_matrix(aMb);
    
    // p4 = (b-d)(g+h)
    Matrix bmd = subtract_matrix(b, d);
    Matrix gMh = add_matrix(g, h);
    Matrix p4 = strassen(bmd, gMh);
    free_matrix(bmd);
    free_matrix(gMh);
    
    // p5 = a(f-h)
    Matrix fmh = subtract_matrix(f, h);
    Matrix p5 = strassen(a, fmh);
    free_matrix(fmh);
    
    // p6 = (c+d)e
    Matrix cMd = add_matrix(c, d);
    Matrix p6 = strassen(cMd, e);
    free_matrix(cMd);
    
    // p7 = (a-c)(e+f)
    Matrix amc = subtract_matrix(a, c);
    Matrix eMf = add_matrix(e, f);
    Matrix p7 = strassen(amc, eMf);
    free_matrix(amc);
    free_matrix(eMf);
    
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
    
    // se unen las 4 submatrices
    Matrix C = join_matrices(C11, C12, C21, C22);
    
    // Liberar memoria de submatrices y productos intermedios
    free_matrix(a); free_matrix(b); free_matrix(c); free_matrix(d);
    free_matrix(e); free_matrix(f); free_matrix(g); free_matrix(h);
    free_matrix(p1); free_matrix(p2); free_matrix(p3); free_matrix(p4);
    free_matrix(p5); free_matrix(p6); free_matrix(p7);
    free_matrix(C11); free_matrix(C12); free_matrix(C21); free_matrix(C22);
    
    // De esta manera T(n) = 7T(n/2) + Theta(n^2) que por teorema maestro es Theta(n^2.81)
    // Esto es mejor que Theta(n^3) de Brute Force, 
    // por lo tanto, el algoritmo de Strassen es más eficiente que el algoritmo de Brute Force.//

    return C;
}

#endif