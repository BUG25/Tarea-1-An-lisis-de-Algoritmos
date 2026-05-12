def brute_force(A,B):
    n1,n2,n3 = A.shape[0], A.shape[1], B.shape[1]  # n1,n2,n3 = filas A, columnas A = filas B, columnas B
    C = np.array([[0]*n3 for i in range(n1)])      # C un array de n1*n3 lleno de ceros
    for i in range(n1):
        for j in range(n2):
            for k in range(n3):
                    C[i][j] += A[i][k]*B[k][j]     # multiplicación normal de matrices, cada espacio ij de C es la suma de multiplicar
                                                   # espacio ik de A por espacio kj de B, por cada k
    return C

