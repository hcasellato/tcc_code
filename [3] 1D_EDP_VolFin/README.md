# Volumes Finitos para problemas unidimensionais

Esses códigos implementam a discretização de volumes finitos para resolver problemas (contínuos e descontínuos) unidimensionais, também resolvem o sistema linear por três métodos de fatoração: de Crout (CF), de Cholesky (Ch) e LDL^t.

- Código `1D_EDP_*_ex1.cpp` implementam as fatorações para o exercício unidimensional contínuo;
- Código `1D_EDP_CF_ex2.cpp` implementa a fatoração de Crout para o exercício unidimensional descontínuo;
- Código `1D_EDP-CF_Pa.cpp` implementa a fatoração de Crout e tenta paralelizar a resolução;
- Pastas `1D_finVolMet_*` são pastas de dados da resolução, para visualização gráfica.