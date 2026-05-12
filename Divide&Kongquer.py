# * RECORDAR: no se pueden usar fokin librerías de multiplicación de matrices! debemos implementarlo nosotro!!!!!!!
# n es múltiplo de 2
# si n=1, se hace la multiplicación directamente

# === ALGORITMO DE STRASSEN ===

# dividimos las matrices en matrices más enanas
def split(matrix):
    n = len(matrix)
    #dividimos en 4 submatrice
    return (matrix[:n//2, :n//2], # 1ra matriz, filas 0 a n//2, columnas 0 a n//2
           matrix[:n//2, n//2:],  # 2da matriz, filas 0 a n//2, columnas n//2 a n
           matrix[n//2:, :n//2],  # 3ra matriz, filas n//2 a n, columnas 0 a n//2
           matrix[n//2:, n//2:])  # 4ta matriz, filas n//2 a n, columnas n//2 a n

