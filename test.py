import matplotlib.pyplot as plt
import numpy as np

# Données pour la heatmap
data = np.random.rand(5, 10)

# Création de la heatmap
plt.imshow(data, cmap='viridis', aspect=0.5)  # Changez "aspect" pour ajuster la forme des cellules
plt.colorbar()
plt.show()