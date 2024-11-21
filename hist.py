import numpy as np
import os
import random
import matplotlib.pyplot as plt

def read_conditions(filename, num_condizioni, neq):
    if not os.path.exists(filename):
        print("Errore nell'apertura del file per lettura!")
        return None

    dimensione_condizione = neq * np.dtype('float64').itemsize
    file_size = os.path.getsize(filename)
    numero_condizioni_tot = file_size // dimensione_condizione

    print("Il numero di condizioni totali è", numero_condizioni_tot)
    if num_condizioni > numero_condizioni_tot:
        print("Errore: il numero di condizioni richiesto eccede il numero di condizioni nel file.")
        return None

    condizioni = np.empty((num_condizioni, neq), dtype='float64')
    with open(filename, 'rb') as inFile:
        selected_indices = random.sample(range(numero_condizioni_tot), num_condizioni)
        for i, indice_casuale in enumerate(selected_indices):
            offset = indice_casuale * dimensione_condizione
            inFile.seek(offset)
            condition_data = np.fromfile(inFile, dtype='float64', count=neq)
            if condition_data.size != neq:
                print("Errore durante la lettura del file!")
                return None
            condizioni[i] = condition_data

    return condizioni

# Input parameters
filename = "condizioni.bin"
num_condizioni = 1000  # Number of conditions to read
N = 10
dim = 1
neq = (N + 2) * 2 * dim + 2

condizioni = read_conditions(filename, num_condizioni, neq)
if condizioni is None:
    exit()

# Analyze position differences
# Estrai x_values: colonne alternate dalla prima, escludendo i termostati
x_values = condizioni[:, :-2:2]

# Estrai p_values: colonne alternate dalla seconda, escludendo i termostati e le particelle fisse
p_values = condizioni[:, 3:-4:2]

chi_values = condizioni[:, -2:]
print(chi_values[:5])

diff_x_all = np.diff(x_values, axis=1).flatten()  # Compute differences and flatten

p_all = p_values.flatten()
# Plot histogram of differences
plt.hist(diff_x_all, bins=20, density=True, alpha=0.7, color='blue')

# Energy and density functions
def potential(x, a, b):
    return (x - 1)**2 / 2 + a * (x - 1)**3 / 3 + b * (x - 1)**4 / 4

def cinetica(x):
    return x**2 / 2

def density_x(x, a, b, beta):
    return np.exp(-beta * potential(x, a, b))

def density_p(x, beta):
    return np.exp(-beta * cinetica(x))


x_vals = np.linspace(min(diff_x_all), max(diff_x_all), 500)
density_vals = density_x(x_vals, 0, 1, 1)
density_vals /= np.trapz(density_vals, x_vals)  # Normalize

p_vals = np.linspace(min(p_all), max(p_all), 500)
density_vals_p = density_p(p_vals-1,1)
density_vals_p /= np.trapz(density_vals_p, p_vals)  # Normalize

# Plot density curve
plt.plot(x_vals, density_vals, color='red', linewidth=2, label='Energy-based Density')

# Add labels and legend
plt.xlabel('Difference between consecutive positions')
plt.ylabel('Probability Density')
plt.title('Histogram of Position Differences with Energy-based Density')
plt.legend()
plt.show()



