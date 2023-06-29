import numpy as np
import matplotlib.pyplot as plt

W = [311924, 305960, 234835, 228871]

plt.plot(W, label = 'Peso da aeronave')
plt.legend()
plt.xlabel('Etapas da missão')
plt.ylabel('Peso (kg)')
plt.show()