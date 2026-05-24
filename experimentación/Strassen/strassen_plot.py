import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Cargar los datos
int_random = pd.read_csv('int_random.csv')
real_random = pd.read_csv('real_random.csv')
identity = pd.read_csv('identity.csv')
constant = pd.read_csv('constant.csv')

# Crear figura
plt.figure(figsize=(10, 6))

# Graficar cada experimento (usando t_mean)
plt.loglog(int_random['n'], int_random['t_mean'], 'o-', label='Enteros aleatorios', linewidth=2, markersize=6)
plt.loglog(real_random['n'], real_random['t_mean'], 's-', label='Reales aleatorios', linewidth=2, markersize=6)
plt.loglog(identity['n'], identity['t_mean'], '^-', label='Identidad', linewidth=2, markersize=6)
plt.loglog(constant['n'], constant['t_mean'], 'd-', label='Constantes', linewidth=2, markersize=6)

# Curva teórica O(n^2.81) ajustada
n_teorico = np.array([1, 33, 65, 97, 129, 161, 193, 225])
# Ajustar factor constante para que coincida aproximadamente con los datos
factor = int_random['t_mean'].iloc[-1] / (n_teorico[-1] ** 2.81)
tiempo_teorico = factor * (n_teorico ** 2.81)
plt.loglog(n_teorico, tiempo_teorico, 'k--', label='Θ(n^{2.81}) teórica', linewidth=2, alpha=0.7)

# Etiquetas y formato
plt.xlabel('Tamaño de matriz (n)', fontsize=12)
plt.ylabel('Tiempo promedio (nanosegundos)', fontsize=12)
plt.title('Algoritmo de Strassen - Comparación de tipos de matrices', fontsize=14)
plt.legend(loc='upper left', fontsize=10)
plt.grid(True, alpha=0.3, linestyle='--')

# Guardar
plt.tight_layout()
plt.savefig('strassen_plot.png', dpi=150)
plt.show()

# Imprimir pendientes para verificar
print("Pendientes aproximadas (log-log):")
for n in [65, 129, 225]:
    pendiente_int = np.log(int_random[int_random['n']==n]['t_mean'].values[0] / int_random[int_random['n']==65]['t_mean'].values[0]) / np.log(n/65)
    print(f"  Enteros n=65→{n}: {pendiente_int:.2f}")