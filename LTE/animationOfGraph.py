import numpy as np
import matplotlib.pyplot as plt
import imageio.v2 as imageio
import os

# Criar pasta para os gráficos, se quiser
output_dir = 'graficos_por_tempo'
os.makedirs(output_dir, exist_ok=True)

# Carregar os dados
data = np.loadtxt("LTE_data/1D_LTE_UW_100.txt", delimiter=';', skiprows=2)

# Separar as colunas
t_all = data[:, 0]       # Coluna 0 -> t
x_all = data[:, 1]       # Coluna 1 -> x
U_all = data[:, 2]       # Coluna 2 -> U(x)
exact_all = data[:, 3]   # Coluna 3 -> exact_solution

# Obter todos os valores únicos de t
tempos_unicos = np.unique(t_all)

for t in tempos_unicos:
    # Filtrar os dados para o tempo atual
    mask = t_all == t
    x = x_all[mask]
    U = U_all[mask]
    exact = exact_all[mask]

    # Plotar
    plt.figure(figsize=(8, 5))
    plt.plot(x, U, label='Solução Numérica U(x)', marker='o', linestyle='-', color='blue')
    plt.plot(x, exact, label='Solução Exata', marker='x', linestyle='--', color='red')
    
    plt.title(f'Comparação U(x) vs Exata - Tempo t = {t:.4f}')
    plt.xlabel('x')
    plt.ylabel('U(x)')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    
    # Salvar a figura
    filename = f't_{t:.4f}.png'
    filepath = os.path.join(output_dir, filename)
    plt.savefig(filepath)
    plt.close()

    print(f'Gráfico salvo: {filepath}')

# Lista os arquivos PNG em ordem de tempo (assumindo que o nome começa com t_...)
image_files = sorted(
    [os.path.join(output_dir, f) for f in os.listdir(output_dir) if f.endswith('.png')],
    key=lambda x: float(x.split('_')[-1].replace('.png', ''))  # Ordena com base no valor de t no nome
)

# Nome do GIF final
output_gif = 'solucao_por_tempo.gif'

# Criar o GIF
with imageio.get_writer(output_gif, mode='I', duration=5, loop=0) as writer:
    for filename in image_files:
        image = imageio.imread(filename)
        writer.append_data(image)
        print(f'Adicionando {filename} ao GIF')

print(f'GIF salvo como: {output_gif}')