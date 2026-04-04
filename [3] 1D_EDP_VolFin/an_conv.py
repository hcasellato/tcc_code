import matplotlib.pyplot as plt

data = [
    (10,     0.01878992369, 0.01878992369, 0.01878992369),
    (100,    0.000251408914, 0.000251408914, 0.000251408914),
    (1000,   0.000002616976, 0.000002616976, 0.000002616976),
    (10000,  0.000000026296, 0.000000026298, 0.000000026297),
    (100000, 0.000000000377, 0.000000000364, 0.000000000178),
]

n = [row[0] for row in data]
ldlt = [row[1] for row in data]
crout = [row[2] for row in data]
cholesky = [row[3] for row in data]

plt.figure(figsize=(8, 5))

plt.plot(n, ldlt, marker='o', label='LDLt')
plt.plot(n, crout, marker='o', label='Crout')
plt.plot(n, cholesky, marker='o', label='Cholesky')

plt.xscale('log')
plt.yscale('log')

plt.xlabel('Nro de Subintervalos (log)')
plt.ylabel('Erro (log)')
plt.title('Análise de Convergência dos Métodos')
plt.legend()
plt.grid(True)

plt.savefig("analise_convergencia.jpeg", dpi=300)