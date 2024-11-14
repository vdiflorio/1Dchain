import numpy as np
import matplotlib.pyplot as plt

# Carica i dati dal file
data = np.loadtxt('ttcf.dat')

# Assegna ogni colonna a una variabile
tempo = data[0:-1:100, 0]
valore1 = data[0:-1:100, 1]
valore2 = data[0:-1:100, 2]

# Crea il grafico
plt.figure()
plt.plot(tempo, valore1, label='Valore 1')
# plt.plot(tempo, valore2, label='Valore 2')

# Etichetta gli assi
plt.xlabel('Tempo')
plt.ylabel('Valore')
plt.title('Grafico dei valori in funzione del tempo')
plt.legend()

# Mostra il grafico
plt.show()
