import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# =====================================
# 1. LEITURA DOS DADOS
# =====================================
df = pd.read_csv("2D_finVolMet_BLUv5_625.txt", sep=";", skiprows=2)

# =====================================
# 2. CRIAR MALHA USANDO PIVOT
# =====================================
X = df.pivot(index="y", columns="x", values="x").values
Y = df.pivot(index="y", columns="x", values="y").values
F = df.pivot(index="y", columns="x", values="f(x,y)").values
U = df.pivot(index="y", columns="x", values="Vx").values
V = df.pivot(index="y", columns="x", values="Vy").values

# =====================================
# 3. PLOTAGEM
# =====================================
plt.figure(figsize=(8, 7))

magnitude = np.sqrt(U**2 + V**2)
U_norm = U / (magnitude + 1e-12)
V_norm = V / (magnitude + 1e-12)

pc1 = plt.pcolormesh(X, Y, F, shading="auto")

plt.quiver(X, Y, U_norm, V_norm, color="white", scale=30)

plt.colorbar(pc1)
plt.tight_layout()
plt.savefig("desenvolvimento_2Deliptica_ex2.png", dpi=300, bbox_inches="tight")
plt.close()
