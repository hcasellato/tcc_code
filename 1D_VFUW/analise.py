import pandas as pd
import matplotlib.pyplot as plt

# Lê o arquivo com separador ';'
df = pd.read_csv("1D_VFUW_dta/1D_VF_100.txt", sep=';', skiprows=3)

# Criação dos gráficos lado a lado
fig, axs = plt.subplots(1, 2, figsize=(12, 5))

# Gráfico 1: x vs pressão (w e y)
axs[0].plot(df['x'], df['w'], marker='o', color='tab:blue', label='w (Aproximado)')
axs[0].plot(df['x'], df['y'], marker='x', color='tab:red', linestyle='--', label='y (Exato)')
axs[0].set_title('x vs vs Pressão (w e y)')
axs[0].set_xlabel('x')
axs[0].set_ylabel('Pressão')
axs[0].legend()
axs[0].grid(True)

# Gráfico 2: x vs velocidade (Vel)
axs[1].plot(df['x'], df['Vel'], marker='s', color='tab:green')
axs[1].set_title('x vs Velocidade')
axs[1].set_xlabel('x')
axs[1].set_ylabel('Velocidade')
axs[1].grid(True)

plt.tight_layout()
plt.savefig("grafico_pressao_vel.png", dpi=300, bbox_inches='tight')
plt.close()