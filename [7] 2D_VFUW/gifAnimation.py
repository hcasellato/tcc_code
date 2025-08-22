import numpy as np
import matplotlib.pyplot as plt
import imageio.v2 as imageio
import os

# =====================================
# 1. LEITURA DOS DADOS
# =====================================
data = np.genfromtxt("2D_VFUW_db/2D_UW_625.txt", delimiter=';', skip_header=2)

# Separar as colunas
Time = data[:, 0]
x = data[:, 1]
y = data[:, 2]
C = data[:, 3]

# Obter dimensões
T  = len(np.unique(Time))
nx = len(np.unique(x))
ny = len(np.unique(y))

print(f"T = {T}, nx = {nx}, ny = {ny}")

# Reorganizar apenas os dados de C
C = C.reshape((T, nx, ny))

# Criar pasta para os frames
output_dir = 'graficos_binarios'
os.makedirs(output_dir, exist_ok=True)

# =====================================
# 4. PLOTAGEM
# =====================================
for t in range(T):
    frame = C[t]

    plt.figure(figsize=(8, 6))
    plt.imshow(frame.T, aspect='auto', origin='lower', cmap='viridis')
    plt.colorbar(label='Concentração')
    plt.title(f'Tempo t = {np.unique(Time)[t]}')
    plt.xlabel('i (M)')
    plt.ylabel('j (N)')
    plt.tight_layout()

    fname = os.path.join(output_dir, f'frame_{t:04d}.png')
    plt.savefig(fname)
    plt.close()
    print(f'Frame salvo: {fname}')

# Criar GIF
output_gif = 'animacao_solucao.gif'
image_files = sorted(
    [os.path.join(output_dir, f) for f in os.listdir(output_dir) if f.endswith('.png')],
    key=lambda x: int(x.split('_')[-1].replace('.png', ''))
)

with imageio.get_writer(output_gif, mode='I', duration=0.5, loop=0) as writer:
    for fname in image_files:
        image = imageio.imread(fname)
        writer.append_data(image)
        print(f'Adicionado ao GIF: {fname}')

print(f'GIF final salvo como: {output_gif}')
