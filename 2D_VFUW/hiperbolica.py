import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.colors import SymLogNorm

# =====================================
# 1. LEITURA DOS DADOS DO CAMPO f(x,y), Vx, Vy
# =====================================
dados = np.genfromtxt("2D_VFUW_db/2D_VF_625.txt", delimiter=';', skip_header=4)

x = dados[:, 0]
y = dados[:, 1]
fxy = dados[:, 2]
vx = dados[:, 3]
vy = dados[:, 4]

# =====================================
# 2. LEITURA DO CAMPO DE PERMEABILIDADE
# =====================================
k_dados = np.genfromtxt("2D_VFUW_db/KField_625.txt", delimiter=';', skip_header=2)
xk = k_dados[:, 0]
yk = k_dados[:, 1]
kxy = k_dados[:, 2]

# =====================================
# 3. REORGANIZAR DADOS EM MALHAS 2D
# =====================================
# Descobrir dimensões da malha automaticamente
x_unique = np.unique(x)
y_unique = np.unique(y)
nx, ny = len(x_unique), len(y_unique)

X, Y = np.meshgrid(x_unique, y_unique)

# Redimensionar os campos
F = fxy.reshape((ny, nx))
K = kxy.reshape((ny, nx))
U = vx.reshape((ny, nx))
V = vy.reshape((ny, nx))

# Normalizar campo vetorial
magnitude = np.sqrt(U**2 + V**2)
U_norm = U / (magnitude + 1e-12)
V_norm = V / (magnitude + 1e-12)

# =====================================
# 4. PLOTAGEM
# =====================================
fig, axs = plt.subplots(1, 2, figsize=(14, 6))

# Heatmap da permeabilidade (log)
pc = axs[0].pcolormesh(X, Y, K, shading='auto', norm=LogNorm())
axs[0].set_title("Permeabilidade (log)")
fig.colorbar(pc, ax=axs[0])

# Heatmap da função f(x,y) com vetores
norma = SymLogNorm(linthresh=1e-3, linscale=1, vmin=np.min(F), vmax=np.max(F))
pc1 = plt.pcolormesh(X, Y, F, shading='auto', norm=norma)
axs[1].quiver(X, Y, U_norm, V_norm, color='white', scale=30)
axs[1].set_title("f(x,y) e campo vetorial normalizado")
fig.colorbar(pc1, ax=axs[1])

plt.tight_layout()
plt.savefig("produto.png", dpi=300, bbox_inches='tight')
plt.close()