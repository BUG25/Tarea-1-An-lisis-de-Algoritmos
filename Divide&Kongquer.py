# * RECORDAR: no se pueden usar fokin librerías de multiplicación de matrices! debemos implementarlo nosotro!!!!!!!
# n es múltiplo de 2
# si n=1, se hace la multiplicación directamente

import BruteForce
import numpy as np

# === ALGORITMO DE STRASSEN ===

# dividimos las matrices en matrices más enanas
def split(matrix):
    n = len(matrix)
    #dividimos en 4 submatrice
    return (matrix[:n//2, :n//2], # 1ra matriz, filas 0 a n//2, columnas 0 a n//2
           matrix[:n//2, n//2:],  # 2da matriz, filas 0 a n//2, columnas n//2 a n
           matrix[n//2:, :n//2],  # 3ra matriz, filas n//2 a n, columnas 0 a n//2
           matrix[n//2:, n//2:])  # 4ta matriz, filas n//2 a n, columnas n//2 a n

def strassen(A,B):
    if len(A) <= 2:
        return BruteForce.brute_force(A,B)  # multiplicamos directamente
    a,b,c,d = split(A)
    e,f,g,h = split(B)     
    # Podríamos calcular con strassen los siguientes productos:
    # ae = strassen(a,e)
    # bg = strassen(b,g)
    # af = strassen(a,f)
    # bh = strassen(b,h)
    # ce = strassen(c,e)
    # dg = strassen(d,g)
    # cf = strassen(c,f)
    # dh = strassen(d,h)
    # y luego C11, C12, C21, C22    
    # y unirlas en una sola matriz C  
    # Tendríamos 8 llamadas a strassen de complejidad T(n/2), y cinco operaciones de complejidad Theta(n^2) 
    # T(n) = 8T(n/2) + Theta(n^2) => T(n) = Theta(n^3) por teorema maestro
    # El problema es que Brute Force tiene la misma complejidad temporal, lo que hace inoca la implementación de algoritmo de Strassen
    # Por tanto, podemos reescribir de manera tal que:
    # 1. ae+bg = ae+ah+de+dh+dg-de-ah-bh+bg+bh-dg-dh = (a+d)(e+h)+d(g-e)-(a+b)h+(b-d)(g+h)
    # 2. af+bh = af-ah+ah+bh = a(f-h)+h(a+b)
    # 3. ce+dg = ce+de+dg-de = e(c+d)+d(g-e)
    # 4. cf+dh = af-ah+ae+ah+de+dh-ce-de-ae-af+ce+cf = a(f-h)+(a+d)(e+h)-e(c+d)-(a-c)(e+f)
    # Con estas expresiones, podemos definir:
    p1 = strassen(a+d,e+h)  # p1 = (a+d)(e+h)
    p2 = strassen(d,g-e)    # p2 = d(g-e)
    p3 = strassen(a+b,h)    # p3 = (a+b)h
    p4 = strassen(b-d,g+h)  # p4 = (b-d)(g+h)
    p5 = strassen(a,f-h)    # p5 = a(f-h)
    p6 = strassen(c+d,e)    # p6 = (c+d)e
    p7 = strassen(a-c,e+f)  # p7 = (a-c)(e+f)
    C11 = p1+p2-p3+p4
    C12 = p5+p3
    C21 = p6+p2
    C22 = p5+p1-p6-p7
    C = np.vstack((np.hstack((C11,C12)), np.hstack((C21,C22))))
    # De esta manera T(n) = 7T(n/2) + Theta(n^2) que por teorema maestro es Theta(n^2.81)
    # Esto es mejor que Theta(n^3) de Brute Force, por lo que el algoritmo de Strassen es más eficiente que el algoritmo de Brute Force
    return C