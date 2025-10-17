import pandas as pd
import matplotlib.pyplot as plt

# Lê o arquivo com separador ';'
df = pd.read_csv("1D_VFUW_dta/1D_VF_100.txt", sep=';', skiprows=3)

# Criação dos gráficos lado a lado
fig, axs = plt.subplots(1, 2, figsize=(12, 5))

# Cores daltônicas seguras
cor_y  = '#0072B2'  # Azul escuro
cor_w  = '#E69F00'  # Amarelo dourado
cor_w2 = '#009E73'

# Gráfico 1: x vs pressão (w e y)
axs[0].plot(df['x'], df['w'], marker='o', label='Solução aproximada $w_i$', color=cor_w, markersize=6)
axs[0].plot(df['x'], df['y'], marker='x', label='Solução exata $y_i$', color=cor_y, markersize=6)
axs[0].set_title('Campo de Pressões')
axs[0].set_xlabel('Posição')
axs[0].set_ylabel('Pressão')
axs[0].legend()
axs[0].grid(True)

# Gráfico 2: x vs velocidade (Vel)
axs[1].plot(df['x'], df['Vel'], marker='s', color=cor_w2, markersize=6)
axs[1].set_title('Campo de Velocidades')
axs[1].set_xlabel('Posição')
axs[1].set_ylabel('Velocidade')
axs[1].grid(True)

plt.tight_layout()
plt.savefig("grafico_pressao_vel.png", dpi=300, bbox_inches='tight')
plt.close()