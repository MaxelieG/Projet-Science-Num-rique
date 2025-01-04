import matplotlib.pyplot as plt

fig, axs = plt.subplots(1, 3, figsize=(10, 4))  # 3 subplots sur une ligne
for ax in axs:
    ax.plot([1, 2, 3], [1, 4, 9])
    ax.set_title("Titre")

fig.tight_layout()  # Ajustement automatique
plt.show()