import numpy as np
import matplotlib.pyplot as plt
import imageio.v2 as imageio
import os

# =====================================
# 1. LEITURA DOS DADOS
# =====================================
data = np.genfromtxt("2D_VFUW_db/debug_VFH.txt", delimiter=';', skip_header=1)

# Separar as colunas
Time = data[:, 0]
x = data[:, 1]
y = data[:, 2]
C = data[:, 3]

x_unique = np.unique(x)
y_unique = np.unique(y)

# Obter dimensões
T  = len(np.unique(Time))
nx = len(x_unique)
ny = len(y_unique)

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

    plt.figure(figsize=(7, 6))
    plt.imshow(
        frame.T,
        aspect='auto',
        origin='lower',
        cmap='viridis',
        extent=[0.0, 1.0, 0.0, 1.0],
        vmin = 0.0,
        vmax = 1.0
    )
    plt.colorbar(label='Concentração')
    plt.title(f'Tempo t = {np.unique(Time)[t]}')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.tight_layout()

    fname = os.path.join(output_dir, f'frame_{t:04d}.png')
    plt.savefig(fname)
    plt.close()

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

print(f'GIF final salvo como: {output_gif}')

# =====================================
# 5. FIGURA COM 6 GRÁFICOS
# =====================================
import matplotlib.gridspec as gridspec

# Escolher 6 índices distribuídos ao longo do tempo
indices = np.linspace(5, T-1, 6, dtype=int)

fig, axes = plt.subplots(2, 3, figsize=(10, 8))

for ax, idx in zip(axes.flat, indices):
    frame = C[idx]

    im = ax.imshow(
        frame.T,
        aspect='auto',
        origin='lower',
        cmap='viridis',
        extent=[0.0, 1.0, 0.0, 1.0],
        vmin=0.0,
        vmax=1.0
    )
    ax.set_title(f"t = {np.unique(Time)[idx]}")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect('equal')

plt.tight_layout(rect=[0.13, 0.1, 1, 0.95])

# Barra de cores embaixo
cbar = fig.colorbar(im, ax=axes.ravel().tolist(), 
                    orientation='horizontal', fraction=0.02, pad=0.02)

plt.savefig("painel_6_graficos.png", dpi=300, bbox_inches="tight")