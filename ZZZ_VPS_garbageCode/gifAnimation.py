import numpy as np
import matplotlib.pyplot as plt
import imageio.v2 as imageio
import os

# === Caminho para o arquivo .bin ===
filename = 'VPS_db/VPS_FVM_625.bin'  # troque o sufixo conforme o DMR

# === Abrir o arquivo binário e ler ===
with open(filename, 'rb') as f:
    # Ler cabeçalho: 3 inteiros (T, M, N)
    dims = np.fromfile(f, dtype=np.int32, count=3)
    T, M, N = dims
    print(f'T = {T}, M = {M}, N = {N}')

    # Ler todo o resto como double
    data = np.fromfile(f, dtype=np.float64)

# Verificar tamanho
expected_size = T * M * N
if data.size != expected_size:
    raise ValueError(f"Tamanho inesperado: esperava {expected_size}, mas li {data.size}")

# Reorganizar dados: (T, M, N)
data = data.reshape((T, M, N))

# === Criar diretório para salvar os gráficos ===
output_dir = 'graficos_binarios'
os.makedirs(output_dir, exist_ok=True)

# === Gerar gráficos por tempo ===
for t in range(T):
    frame = data[t]  # shape (M, N)

    plt.figure(figsize=(8, 6))
    plt.imshow(frame.T, aspect='auto', origin='lower', cmap='viridis', vmin=0, vmax=1)
    plt.colorbar(label='Valor')
    plt.title(f'Tempo t = {t}')
    plt.xlabel('i (M)')
    plt.ylabel('j (N)')
    plt.tight_layout()

    fname = os.path.join(output_dir, f'frame_{t:04d}.png')
    plt.savefig(fname)
    plt.close()
    print(f'Frame salvo: {fname}')

# === Gerar GIF ===
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
