import numpy as np
import matplotlib.pyplot as plt

# =====================================
# 1. LEITURA DOS DADOS
# =====================================
dados = np.genfromtxt("2D_VFUW_db/debug_VFE.txt", delimiter=';', skip_header=1)

x = dados[:, 0]
y = dados[:, 1]
pressao = dados[:, 2]   # f(x,y) = pressão
vx = dados[:, 3]
vy = dados[:, 4]

# =====================================
# 2. REORGANIZAR DADOS EM MALHAS 2D
# =====================================
x_unique = np.unique(x)
y_unique = np.unique(y)
nx, ny = len(x_unique), len(y_unique)

X, Y = np.meshgrid(x_unique, y_unique)

P = pressao.reshape((ny, nx))
U = vx.reshape((ny, nx))
V = vy.reshape((ny, nx))

# Normalizar campo vetorial
magnitude = np.sqrt(U**2 + V**2)
U_norm = U / (magnitude + 1e-12)
V_norm = V / (magnitude + 1e-12)

# =====================================
# 3. PLOTAGEM
# =====================================
plt.figure(figsize=(7,6))

pc = plt.pcolormesh(X, Y, P, cmap="viridis")
plt.colorbar(pc, label="Pressão")

plt.quiver(X, Y, U_norm, V_norm, color='white', scale=50)
plt.title("Pressão e Campo de Velocidade")
plt.xlabel("x")
plt.ylabel("y")

plt.tight_layout()
plt.savefig("pressao_velocidade.png", dpi=300, bbox_inches='tight')
plt.close()
