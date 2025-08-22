import pandas as pd
import matplotlib.pyplot as plt
import imageio.v2 as imageio
import os

# Leitura dos dados
df = pd.read_csv("1D_VFUW_dta/1D_UW_100.txt", sep=';', skiprows=2)

# Criar pasta para armazenar os frames temporários
folder = "frames_temp"
os.makedirs(folder, exist_ok=True)

# Obter todos os tempos únicos ordenados
tempos = sorted(df['t'].unique())

# Gerar um gráfico para cada tempo
filenames = []
for i, t_val in enumerate(tempos):
    subset = df[df['t'] == t_val]
    x = subset['x']
    Co = subset['Co']

    # Separar os dados em dois grupos: Co ≤ 1 e Co > 1
    mask_normal = Co <= 1.0
    mask_high = Co > 1.0

    plt.figure(figsize=(6, 4))

    # Pontos com Co ≤ 1.0: azul
    plt.scatter(x[mask_normal], Co[mask_normal], color='blue', label='Co ≤ 1')

    # Pontos com Co > 1.0: vermelho
    plt.scatter(x[mask_high], Co[mask_high], color='red', label='Co > 1')

    plt.xlabel("Posição x")
    plt.ylabel("Concentração Co")
    plt.title(f"Tempo t = {t_val}")
    plt.grid(True)
    plt.xlim(0, 1.0)
    plt.ylim(0, 1.3)
    filename = f"{folder}/frame_{i:03d}.png"
    plt.savefig(filename)
    filenames.append(filename)
    plt.close()

# Criar o GIF
with imageio.get_writer("concentracao_x_posicao.gif", mode='I', duration=0.3, loop = 0) as writer:
    for filename in filenames:
        image = imageio.imread(filename)
        writer.append_data(image)

# (Opcional) Remover os frames temporários
for filename in filenames:
    os.remove(filename)
os.rmdir(folder)

print("GIF criado com sucesso: concentracao_x_posicao.gif")
