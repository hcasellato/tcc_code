import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# =====================================
# 1. LEITURA DO CAMPO DE PERMEABILIDADE
# =====================================
k_dados = np.genfromtxt("2D_VFUW_db/KField_625.txt", delimiter=';', skip_header=2)
xk = k_dados[:, 0]
yk = k_dados[:, 1]
kxy = k_dados[:, 2]

# =====================================
# 2. REORGANIZAR DADOS EM MALHA 2D
# =====================================
x_unique = np.unique(xk)
y_unique = np.unique(yk)
nx, ny = len(x_unique), len(y_unique)

X, Y = np.meshgrid(x_unique, y_unique)
K = kxy.reshape((ny, nx))

# =====================================
# 3. PLOTAGEM
# =====================================
plt.figure(figsize=(7, 6))
plt.title("Campo de permeabilidades (log)")
plt.xlabel("x", rotation=0)
plt.ylabel("y", rotation=0)

pcm = plt.pcolormesh(X, Y, K, shading='auto', norm=LogNorm(), cmap='viridis')
cbar = plt.colorbar(pcm)
cbar.set_label("K", rotation=0)

plt.tight_layout()

plt.savefig("perm_hetero.png", dpi=300, bbox_inches='tight')
plt.close()
