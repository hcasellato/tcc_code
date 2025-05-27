import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Load the data
file_path = '2D_finVolMet_BLU/debugMatrix.txt'  # Replace with your actual filename
df = pd.read_csv(file_path, sep=';')
df.columns = df.columns.str.strip().str.lower()

# Verifica colunas
assert {'x', 'y', 'm', 'l', 'u'}.issubset(df.columns), "Colunas incorretas"

# Ordena
x_vals = np.sort(df['x'].unique())
y_vals = np.sort(df['y'].unique())[::-1]

# Matrizes
M = df.pivot(index='y', columns='x', values='m').reindex(index=y_vals, columns=x_vals)
L = df.pivot(index='y', columns='x', values='l').reindex(index=y_vals, columns=x_vals)
U = df.pivot(index='y', columns='x', values='u').reindex(index=y_vals, columns=x_vals)

# Escala comum de cores
vmin = min(M.min().min(), L.min().min(), U.min().min())
vmax = max(M.max().max(), L.max().max(), U.max().max())

# Função para plotar com máscara para zeros
def plot_heatmap(matrix, title, ax):
    mask = matrix == 0
    sns.heatmap(matrix, mask=mask, cmap="viridis", vmin=vmin, vmax=vmax,
                cbar=True, ax=ax, square=True)
    ax.set_title(title)
    ax.invert_yaxis()

# Plotagem
fig, axes = plt.subplots(1, 3, figsize=(18, 6))

plot_heatmap(M, 'Matriz M (≠ 0 destacados)', axes[0])
plot_heatmap(L, 'Matriz L (≠ 0 destacados)', axes[1])
plot_heatmap(U, 'Matriz U (≠ 0 destacados)', axes[2])

plt.tight_layout()
fig.savefig('matrizes_destacadas.png', dpi=300)