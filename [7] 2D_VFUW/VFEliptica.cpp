/* ================================ | Código Exercício Chemetov-Henrique | ===================================== */
//
// Implementação de Block LU factorization para matrizes em bloco tridiagonais para os problemas bidimensionais
// do livro de Volumes Finitos (Fabricio). Consultar também ANALYSIS OF NUMERICAL METHODS, Isaacson e Keller.
// 
// Essa versão é a refatoração de 2D_EDP_BLU_v4, tentando otimizar espaço e performance, além de fazer o código
// mais legível.
//
// Problemas da forma: 
//
// NABLA . (- K NABLA p) = q   em OMEGA
//                    p  = p_b em OMEGA_p
//    (- K NABLA p) . n  = u_b em OMEGA_u
// 
// com pressão relacionada à vel. de Darcy
// 
// u = (- K NABLA p)
//
/* ============================================================================================================= */

// Bibliotecas
#include <iostream>
#include <unistd.h>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <string>
#include <cmath>
#include <ctime>

#include <chrono>

using namespace std;


double q(double x, double y) {
  int test   = 1;
  if(x == 0.02*test && y == 0.02*test)
    return  1.0;
  else if(x == 1 - 0.02*test && y == 1 - 0.02*test)
    return -1.0;
  else
    return  0.0;
}

void Copy_A_to_Atil(double **A, double **Atil, int block_index, int TB){
  int start = (block_index - 1) * TB; // DMR - TB = 600

  Atil[1][1] = A[2][start + 1]; // 601
  Atil[1][2] = A[3][start + 1];

  for(int i = 2; i < TB; i++){
    Atil[i][i-1] = A[1][start + i];
    Atil[i][i]   = A[2][start + i];
    Atil[i][i+1] = A[3][start + i];
  }

  Atil[TB][TB-1] = A[1][start + TB]; // 625
  Atil[TB][TB]   = A[2][start + TB];

}

void LU_solve(double **D, double *X, double *F, double **L, double **U, int TB){
  double sum;
  double* z = new double[TB+1];

  // Passo 1
  L[1][1] = D[1][1];

  for (int j = 1; j <= TB; j++) {
    // Calculate column j of L
      for (int i = j; i <= TB; i++) {
        sum = 0.0;
        for (int k = 1; k < j; k++) {
          sum += L[i][k] * U[k][j];
        }
      L[i][j] = D[i][j] - sum;
    }
    
    // Calculate row j of U
    for (int i = j + 1; i <= TB; i++) {
      sum = 0.0;
      for (int k = 1; k < j; k++) {
        sum += L[j][k] * U[k][i];
      }
      U[j][i] = (D[j][i] - sum) / L[j][j];
    }
  }

  // Passo 2 >> Solucionar Lz = F
  z[1] = F[1] / L[1][1];

  // Passo 7.2
  for(int i = 2; i <= TB; i++){
    sum = 0.0;
    for(int j = 1; j <= i - 1; j++)
      sum += L[i][j] * z[j];
    z[i] = (F[i] - sum)/L[i][i];
  }

  // Passo 3 >> Substituição regressiva
  X[TB] = z[TB];

  // Passo 3.1
  for(int i = TB - 1; i >= 1; i--){
    sum = 0.0;
    for(int j = i+1; j <= TB; j++)
      sum += U[i][j] * X[j];
    X[i] = (z[i] - sum);
  }

  // Cleanup memory
  delete[] z;
}

double** VFE(int M, int N, int TB, string name_file, int type_K){
  // ================= | Variáveis!
  int DMR;  // Dimensão da Matriz Resultante (DMR)

  double Al, Be, Ga, De;   // (x,y) \in [Al.Be] X [Ga,De]
  double hx, hy, hx2, hy2; // Pulo entre x_i e x_{i+1}
  
  // =================================================================
  // COLOQUE AQUI OS VALORES DE INTERESSE
  
  Al = Ga = 0.0;
  Be = De = 1.0;

  DMR = M * N;

  // =================================================================
  // DIAGONAIS
  double** diag = new double*[5];
  for(int i = 0; i < 5; i++)
    diag[i] = new double[DMR+1];

  int k; // aritmética k = i + (j - 1)M
  double x, y;
  
  // RESPOSTA!
  double w[DMR + 1];
  
  // Específico
  double Khx = 0.0;
  double Khy = 0.0;

  double** K_cell = new double*[M+1];
  for(int i = 0; i <= M; i++)
    K_cell[i] = new double[M+1];
  
  // =================================================================
  // MATRIZES PARA O BLU
  double*** Atil = new double**[TB+1];
  double*** G    = new double**[TB+1];
  double*** B    = new double**[TB+1];

  double**  C    = new double*[TB+1];
  double**  Z    = new double*[TB+1];

  double**  b_d  = new double*[TB+1]; // bloco de d
  double**  b_w  = new double*[TB+1]; // bloco de w

  double*   b_i  = new double[TB+1]; // vetor intermediario

  // Alocando número de blocos
  for(int i = 0; i <= TB; i++){
    Atil[i] = new double*[TB+1];
    G[i]    = new double*[TB+1];
    B[i]    = new double*[TB+1];

    // Alocando linhas
    C[i]    = new double[TB+1];
    Z[i]    = new double[TB+1];
    
    b_d[i]  = new double[TB+1];
    b_w[i]  = new double[TB+1];

    // Alocando colunas
    for(int j = 0; j <= TB; j++){
      Atil[i][j] = new double[TB+1];
      G[i][j]    = new double[TB+1];
      B[i][j]    = new double[TB+1];
      
      // Inicializando colunas
      C[i][j]    = 0.0;
      Z[i][j]    = 0.0;
      
      b_d[i][j]  = 0.0;
      b_w[i][j]  = 0.0;
      
      // Inicializando colunas
      for(int k = 0; k <= TB; k++){
        Atil[i][j][k] = 0.0; 
        G[i][j][k]    = 0.0; 
        B[i][j][k]    = 0.0; 
      }
    }
  }

  // ============================| Começo 
  // Passo 1
  hx = (Be - Al)/M;
  hy = (De - Ga)/N;

  hx2 = hx * hx;
  hy2 = hy * hy;

  for(int i = 1; i <= M; i++){
    for(int j = 1; j <= N; j++){
      k = i + (j - 1)*M;

      double x = Al + (i - 0.5) * hx;
      double y = Ga + (j - 0.5) * hy;

      K_cell[i][j]  = K_function(x, y, type_K);
      b_d[i][j]     = q(x, y);
    }
  }

  // Passo 2 >> Montagem da Matriz
  for(int i = 1; i <= M; i++){
    for(int j = 1; j <= N; j++){
      k = i + (j - 1)*M;
      
      double left, right, bottom, top;
      
      // Calculate interface permeabilities using harmonic mean
      left   = (i > 1) ? K_half(K_cell[i-1][j], K_cell[i][j]) / hx2 : 0.0;
      right  = (i < M) ? K_half(K_cell[i+1][j], K_cell[i][j]) / hx2 : 0.0;
      bottom = (j > 1) ? K_half(K_cell[i][j-1], K_cell[i][j]) / hy2 : 0.0;
      top    = (j < N) ? K_half(K_cell[i][j+1], K_cell[i][j]) / hy2 : 0.0;

      diag[0][k] = -bottom;    // bottom
      diag[1][k] = -left;      // left
      diag[3][k] = -right;     // right
      diag[4][k] = -top;       // top
      diag[2][k] = left + right + bottom + top;  // center
    }
  }

  // Passo 2.4 >> Introduzindo solução de Dirichlet local em (1,1)
  diag[2][1]              = 1.0; // Diagonal principal igual a 1
  diag[4][1] = diag[3][1] = 0.0; // Resto zerado
  b_d[1][1]   =  1.0;
  
  diag[2][DMR]                = 1.0; // Diagonal principal igual a 1
  diag[0][DMR] = diag[1][DMR] = 0.0; // Resto zerado
  b_d[TB][TB] = -b_d[1][1];
  

  // ============================| Block LU Decomposition
  auto start = std::chrono::high_resolution_clock::now();
  
  // Passo 0 >> Fazer matrizes L e U
  double** L = new double*[TB + 1];
  double** U = new double*[TB + 1];
  for(int i = 0; i <= TB; i++){
    L[i] = new double[TB+1];
    U[i] = new double[TB+1];
    for(int j = 1; j <= TB; j++){
      L[i][j] = 0.0;
      U[i][j] = 0.0;
    }
    U[i][i] = 1.0; // diagonal = 1
  }

  // Passo 3.1.1 >> Ã1 = A1
  Copy_A_to_Atil(diag, Atil[1], 1, TB);

  // Passo 3.1.2 >> Ã1 * G1 = C1 => G1
  for(int i = 1; i <= TB; i++){
    // Pré: fazer C1 (vou ter q fazer isso sempre)
    C[i][i] = diag[4][i];

    LU_solve(Atil[1], G[1][i], C[i], L, U, TB);
  }
  auto end = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> elapsed = end - start;
  std::cout << "Factorization time: " << elapsed.count() << " seconds" << std::endl;

  // Passo 4.1 >> Ã1 * z1 = d1 => z1
  LU_solve(Atil[1], Z[1], b_d[1], L, U, TB);

  // Passo 3.2.1 >> Ãi = Ai - Bi * G[i-1], i = 2, 3, ..., n-1
  for(int b = 2; b < TB; b++){
    int bli = (b-1)*TB;

    // Primeiro copio os valores de diag para Atil no bloco
    Copy_A_to_Atil(diag, Atil[b], b, TB);

    for(int i = 1; i <= TB; i++){
      for(int j = 1; j <= TB; j++){
        // Depois subtraio
        Atil[b][i][j] -= diag[0][bli+i] * G[b-1][i][j];
      }
    }

    // Passo 3.3 >> Ãi * Gi = Ci => Gi,    i = 2, 3, ..., n-1
    for(int i = 1; i <= TB; i++){
      // Pré: fazer Ci
      C[i][i] = diag[4][bli+i];

      LU_solve(Atil[b], G[b][i], C[i], L, U, TB);
    }

    // Passo 4.2 >> Ãi * zi = di - Bi * z[i-1] => zi, i = 2, 3, ..., n-1
    // Fazendo di - Bi * z[i-1]
    for(int j = 1; j <= TB; j++)
      b_i[j] = b_d[b][j] - diag[0][bli+j] * Z[b-1][j];

    // Resolvendo Ãi * zi = bi = di - Bi * z[i-1]
    LU_solve(Atil[b], Z[b], b_i, L, U, TB);
  }

  // Passo 3.2.2 >> Ãn = An - Bn * G[n-1], i = n
  int bli = (TB-1)*TB;

  Copy_A_to_Atil(diag, Atil[TB], TB, TB);

  for(int i = 1; i <= TB; i++){
    for(int j = 1; j <= TB; j++){
      Atil[TB][i][j] -= diag[0][bli+i] * G[TB-1][i][j];
    }
  }

  // Passo 4.2 >> Ãi * zi = di - Bi * z[i-1] => zi, i = n
  // Fazendo di - Bi * z[i-1]
  for(int j = 1; j <= TB; j++)
    b_i[j] = b_d[N][j] - diag[0][bli+j] * Z[N-1][j];

  // Resolvendo Ãi * zi = bi = di - Bi * z[i-1]
  LU_solve(Atil[N], Z[N], b_i, L, U, TB);

  // Passo 5.1 >> xn = zn
  for(int j = 1; j <= TB; j++)
    b_w[TB][j] = Z[TB][j];

  // Passo 5.2 >> xi = zi - Gi * x[i+1],     i = n-1, n-2, ..., 1
  for(int b = TB-1; b >= 1; b--){
    for(int i = 1; i <= TB; i++){
      // Calcular Gi * x[i+1]
      double sum = 0.0;
      for(int j = 1; j <= TB; j++)
        sum += G[b][i][j] * b_w[b+1][j];

      // xi = zi - Gi * x[i+1]
      b_w[b][i] = Z[b][i] - sum;
    }
  }

  // Passo 5.3 >> Transformar block_w to w:
  for(int i = 1; i <= TB; i++){
    int bli = (i-1)*TB;

    for(int j = 1; j <= TB; j++)
      w[bli+j] = b_w[i][j];
  }

  // Passo 5.4 >> Fazer uma matriz de pressao
  double** p = new double*[M+1];
  for(int i = 0; i <= M; i++)
    p[i] = new double[N+1];

  for(int i = 1; i <= M; i++){
    for(int j = 1; j<= N; j++){
      k = i + (j - 1)*M;
      p[i][j] = w[k];
    }
  }

  // Passo 6 >> Campo de Velocidades
  double** campoVel = new double*[2];
  for(int i = 0; i < 2; i++)
    campoVel[i] = new double[DMR+1];

  for(int i = 1; i <= M; i++){
    for(int j = 1; j <= N; j++){
      k = i + (j - 1)*M;

      // Inicializando com 0
      campoVel[0][k] = campoVel[1][k] = 0.0;

      // Passo 6.1
      Khx = (i > 1) ? K_half(K_cell[i][j], K_cell[i-1][j]) : 0.0;
      campoVel[0][k] += -((Khx)* (w[k] - w[k-1])) / hx ;
      
      Khx = (i < M) ? K_half(K_cell[i][j], K_cell[i+1][j]) : 0.0;
      campoVel[0][k] += -((Khx)* (w[k+1] - w[k])) / hx;

      // Passo 6.2
      Khy = (j > 1) ? K_half(K_cell[i][j], K_cell[i][j-1]) : 0.0;
      campoVel[1][k] += -((Khy) * (w[k] - w[k-M])) / hy;

      Khy = (j < N) ? K_half(K_cell[i][j], K_cell[i][j+1]) : 0.0;
      campoVel[1][k] += -((Khy) * (w[k+M] - w[k])) / hy;
      
      campoVel[0][k] /= 2.0;
      campoVel[1][k] /= 2.0;
    }
  }

  // ============================| Print
  // Impressão dos Valores
  // Por favor, crie a pasta antes!
  ofstream file;
  file.open(name_file);

  file << fixed << setprecision(12);
  file << "x;y;f(x,y);Vx;Vy" << endl;

  file << fixed << setprecision(12);
  for(int i = 1; i <= M; i++){
    for(int j = 1; j <= N; j++){
      k = i + (j - 1)*M;
      file << Al + (i - 0.5) * hx  << ";";
      file << Ga + (j - 0.5) * hy  << ";";
      file << p[i][j]              << ";";
      file << campoVel[1][k]       << ";";
      file << campoVel[0][k]       << endl;
    }
  }

  file.close();
  cout << "VFE Completa!" << endl;

  // ============================| Cleanup

  // Cleanup diag array
  for(int i = 0; i < 5; i++)
      delete[] diag[i];
  delete[] diag;

  // Cleanup Atil, G, and B arrays
  for(int i = 0; i <= TB; i++) {

    // Cleanup Atil[i], G[i], B[i]
    if(Atil[i] != nullptr) {
      for(int j = 0; j <= TB; j++)
        delete[] Atil[i][j];
      delete[] Atil[i];
    }
    if(G[i] != nullptr) {
      for(int j = 0; j <= TB; j++)
        delete[] G[i][j];
      delete[] G[i];
    }
    if(B[i] != nullptr) {
      for(int j = 0; j <= TB; j++)
        delete[] B[i][j];
      delete[] B[i];
    }

    // Cleanup C[i], Z[i], b_d[i], b_w[i]
    delete[] C[i];
    delete[] Z[i];
    delete[] b_d[i];
    delete[] b_w[i];
    delete[] L[i];
    delete[] U[i];
    delete[] K_cell[i];
  }

  for (int i = 0; i < 2; i++)
    delete[] campoVel[i];
  delete[] campoVel;

  // Delete the top-level arrays
  delete[] Atil;
  delete[] G;
  delete[] B;
  delete[] C;
  delete[] Z;
  delete[] b_d;
  delete[] b_w;
  delete[] b_i;
  delete[] L;
  delete[] U;
  delete[] K_cell;

  return p;
}

// g++ VFEliptica.cpp -lm -o VFEliptica && ./VFEliptica