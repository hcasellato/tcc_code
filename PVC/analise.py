import matplotlib.pyplot as plt
import numpy as np

# Lê o arquivo
filename = "dta_PVC_TDMA.txt"
data = np.genfromtxt(filename, delimiter=';', skip_header=1)

# Separa os dados
x_i = data[:, 0]
w_i = data[:, 1]
y_i = data[:, 2]

# Cores daltônicas seguras
cor_y = '#0072B2'  # Azul escuro
cor_w = '#E69F00'  # Amarelo dourado

# Gráfico de comparação
plt.figure(figsize=(8, 5))

plt.plot(x_i, y_i, 'x-', label='Solução exata $y_i$', color=cor_y, markersize=6)
plt.plot(x_i, w_i, 'o--', label='Solução aproximada $w_i$', color=cor_w, markersize=6)

plt.xlabel('$x_i$')
plt.ylabel('$y_i$, $w_i$')
plt.title('Comparação entre $y_i$ (exata) e $w_i$ (aproximada)')
plt.legend()
plt.grid(True)

plt.tight_layout()

# Salva a figura
plt.savefig("desenvolvimento_PVC_comparacao.jpeg", dpi=300)