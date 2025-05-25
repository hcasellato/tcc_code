/* ================================ | Código Exercício Chemetov-Henrique | ===================================== */
//
// Tenta a implementação de Block LU factorization para matrizes em bloco tridiagonais para os problemas bidimen-
// sionais do livro de Volumes Finitos. Consultar ANALYSIS OF NUMERICAL METHODS, Isaacson e Keller.
// 
// Essa é minha segunda tentativa. Que os deuses tenham mais misericórdia desta vez.
// 
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
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <string>
#include <cmath>
#include <ctime>

using namespace std;
const double PI  = 3.141592653589793238463;    //value of pi
const double PI2 = 9.869604401089358618834;    //value of pi^2

// =================================================================
// PROBLEMA:
//
// NABLA . u          = -8 PI^2 cos(2PI x)cos(2PI y) em [0,1] X [0,1]
// (- K NABLA p)      = u                            em OMEGA_p
// (- K NABLA p) . n  = 0                            em borda OMEGA
// 
// com permeabilidade absoluta constante
// 
// K = 1
// 
// e solução exata
// 
// p(x,y) = cos(2PI x)cos(2PI y)
//
// Esse problema tem condição de Neumann HOMOGÊNEA, dado que
// u . n = 0 -> fluxo nulo!
// 
// Para consultar implementação heterogênea, consultar livro-txt
// Métodos de Volumes Finitos (pg. 40). 
//
// =================================================================
double K_function(double x) {
  return 1.0;
}

int test_x = 1;
int test_y = 1;

double q(double x, double y) {
  return (4.0 * test_x + 4.0 * test_y) * PI2 * cos(2.0 * PI * x * test_x) * cos(2.0 * PI * y * test_y);
}

// Funciona para ambos 'i' e 'j'
double K_half(double K_1, double K_2) {
  return 2 * (K_1 * K_2) / (K_1 + K_2);
}

// Função exata, se houver
double exact_solution(double x, double y) {
  return cos(2.0 * PI * x * test_x) * cos(2.0 * PI * y * test_y);
}

/*
  0..n = 0 até n INCLUSO!!!

  A = diag[0..4][0..DMR]
  Atil[0..25][0..3][0..25]
*/
void Copy_A_to_Atil(double **A, double **Atil, int block_index){
  int TB = 25;
  int start = (block_index - 1) * TB;

  for(int i = 1; i <= TB; i++){
    Atil[3][i] = A[3][start + i]; // Upper diag
    Atil[2][i] = A[2][start + i]; // Center diag
    Atil[1][i] = A[1][start + i]; // Lower diag
  }
}

// Decompõe e resolve matrizes tridiagonais por decomposição LU
// Dx = f => L(Ux) = f => Lz = f => Ux = z
/*
  INPUT:
    Matriz  D[0..TB]
    Vetor   F[0..TB]
    Vetor   X[0..TB]
*/
void LU_solve(double **D, double* X, double* F){
  int TB = 25;

  double* l = new double[TB + 1];
  double* u = new double[TB + 1];
  double* z = new double[TB + 1];

  // Passo 1.1 >> Decomposição LU
  l[1] = D[2][1];
  u[1] = D[3][1] / l[1];
  z[1] = F[1] / l[1];
  
  // Passo 1.2
  for(int i = 2; i < TB; i++)
  {
    l[i] = D[2][i] - (D[1][i] * u[i-1]);
    u[i] = D[3][i] / l[i];
    z[i] = (F[i] - (D[1][i] * z[i-1])) / l[i];
  }

  // Passo 1.3
  l[TB] = D[2][TB] - (D[1][TB] * u[TB-1]);
  z[TB] = (F[TB]   - (D[1][TB] * z[TB-1])) / l[TB];
  
  // Passo 2.1 >> Passos para 'backward substitution'
  X[TB] = z[TB];

  // Passo 2.2
  for(int i = TB - 1; i >= 1; i--)
    X[i] = z[i] - u[i] * X[i+1];

  delete[] l;
  delete[] u;
  delete[] z;
}

void FinVol()
{
  // ================= | Variáveis!
  int M, N; // Dimensões M,N da matriz e
  int DMR;  // Dimensão da Matriz Resultante (DMR)
  int TB;   // Tamanho do bloco

  double A, Bi, Ci, D;       // (x,y) \in [A.B] X [C,D]
  double alpha, beta;      // Soluções de contorno
  double hx, hy, hx2, hy2; // Pulo entre x_i e x_{i+1}
  
  int debug_const = 1;     // Altere para 1 caso queira ver o tempo de
                           // execução de cada passo
 
  clock_t time_req, inter_time; // Para medir o tempo de execução!
  
  // =================================================================
  // COLOQUE AQUI OS VALORES DE INTERESSE
  
  A = Ci = 0.0;
  Bi = D = 1.0;
  
  M   = 25;
  N   = 25;
  TB  = 25;

  DMR = M * N;

  // =================================================================
  // DIAGONAIS
  /*
    indices para matriz A3x3
  1: p11 p21 p31 
  2: p12 p22 p32 
  3: p13 p23 p33 

  indices da matriz diag A9x9
      _____     distância M (obs. x entre 1...M)
     |     |                            k = x + (y - 1).M
  1: 2 3 . 4 . . . . .   p11       d1   k = 1 + (1 - 1).3 = 1 + 0
  2: 1 2 3 . 4 . . . .   p21       d2   k = 2 + (1 - 1).3 = 2 + 0
  3: . 1 2 . . 4 . . .   p31       d3   k = 3 + (1 - 1).3 = 3 + 0
  4: 0 . . 2 3 . 4 . .   p12       d4   k = 1 + (2 - 1).3 = 1 + 3
  5: . 0 . 1 2 3 . 4 .   p22   =   d5   k = 2 + (2 - 1).3 = 2 + 3
  6: . . 0 . 1 2 . . 4   p32       d6   k = 3 + (2 - 1).3 = 3 + 3
  7: . . . 0 . . 2 3 .   p13       d7   k = 1 + (3 - 1).3 = 1 + 6
  8: . . . . 0 . 1 2 3   p23       d8   k = 2 + (3 - 1).3 = 2 + 6
  9: . . . . . 0 . 1 2   p33       d9   k = 3 + (3 - 1).3 = 3 + 6
               |     |
               |_____| distância M (obs. x entre 1...M)

  lembrando que é de baixo para cima, da esquerda para a direita
  */

  double** diag = new double*[5];
  for(int i = 0; i < 5; i++)
    diag[i] = new double[DMR+1];

  int k; // aritmética k = i + (j - 1)M
  double x, y;
  double d[DMR+1];
  double valor_real[DMR+1];   // Vetor de valores reais da solução
  
  // RESPOSTA!
  double w[DMR + 1];
  
  // Específico
  double Kx[DMR + 2];
  double Ky[DMR + 2];

  // ============================| Começo 
  time_req = inter_time = clock();

  // Passo 1
  hx = (Bi - A)/M;
  hy = (D - Ci)/N;

  hx2 = hx * hx;
  hy2 = hy * hy;
  
  // Passo 1.1 >> Vetor de Ks
  Kx[0]     = K_function(A);
  Kx[DMR+1] = K_function(Bi);

  Ky[0]     = K_function(Ci);
  Ky[DMR+1] = K_function(D);

  for(int i = 1; i <= M; i++){
    for(int j = 1; j <= N; j++){
      k = i + (j - 1)*M;

      Kx[k] = K_function(A + (i - 0.5) * hx);
      Ky[k] = K_function(Ci + (j - 0.5) * hy);

      // Passo 1.2 >> Vetor 'd' e das soluções exatas 
      d[k]          = q(A + (i - 0.5) * hx, Ci + (j - 0.5) * hy);
      valor_real[k] = exact_solution(A + (i - 0.5) * hx, Ci + (j - 0.5) * hy);
    }
  }

  if(debug_const == 1)
    cout << "Passo 1: " << (double)(clock() - inter_time)/CLOCKS_PER_SEC << endl;
  inter_time = clock();

  // Passo 2 >> Montagem da Matriz
  for(int i = 1; i <= M; i++){
    for (int j = 1; j <= N; j++){
      k = i + (j - 1)*M;

      // Passo 2.1 >> Diagonais exteriores (y-direction)
      diag[0][k] = (j > 1) ? -K_half(Ky[k], Ky[k - M]) / hy2 : 0.0;
      diag[4][k] = (j < N) ? -K_half(Ky[k], Ky[k + M]) / hy2 : 0.0; // Fixed M to N

      // Passo 2.2 >> Diagonais internas (x-direction)
      diag[1][k] = (i > 1) ? -K_half(Kx[k], Kx[k - 1]) / hx2 : 0.0;
      diag[3][k] = (i < M) ? -K_half(Kx[k], Kx[k + 1]) / hx2 : 0.0;

      // Passo 2.3 >> Diagonal Central
      diag[2][k] = - diag[0][k] - diag[1][k] - diag[3][k] - diag[4][k];
    }
  }

  if(debug_const == 1)
    cout << "Passo 2: " << (double)(clock() - inter_time)/CLOCKS_PER_SEC << endl;
  inter_time = clock();

  // ============================| Block LU Decomposition
  /*
    [ 2 3 .   4 . .   . . . ]   [      |      |      ]
    [ 1 2 3   . 4 .   . . . ]   [  A1  |  C1  |      ]
    [ . 1 2   . . 4   . . . ]   [      |      |      ]

    [ 0 . .   2 3 .   4 . . ]   [      |      |      ]
    [ . 0 .   1 2 3   . 4 . ] = [  B2  |  A2  |  C2  ] = A
    [ . . 0   . 1 2   . . 4 ]   [      |      |      ]
             
    [ . . .   0 . .   2 3 . ]   [      |      |      ]
    [ . . .   . 0 .   1 2 3 ]   [      |  B3  |  A3  ]
    [ . . .   . . 0   . 1 2 ]   [      |      |      ]

    Ax = f => L(Ux) = f => Ly = f => Ux = y,

    [ Ã1       ] [ I1 G1    ]
    [ B1 Ã2    ] [    I2 G2 ] = LU
    [    B2 Ã3 ] [       I3 ]

  Procedimento:
    
    Ã1 = A1                                         3.1.1
    
    G1 = A1` * C1                 (.` é a inversa)
    Ã1 * G1 = C1 => G1                              3.1.2
    
    Ãi = Ai - Bi * G[i-1], i = 2, 3, ..., n-1       3.2.1
    Gi = Ãi` * Ci
    Ãi * Gi = Ci => Gi,    i = 2, 3, ..., n-1       3.3
    
    Ãn = An - Bn * G[n-1], i = n                    3.2.2

    então,

    z1 = Ã1` * d1
    Ã1 * z1 = d1 => z1                              4.1
    
    zi = Ãi` * (di - Bi * z[i-1]), i = 2, 3, ..., n
    Ãi * zi = di - Bi * z[i-1] => zi                4.2
        
    xn = zn                                         5.1
    xi = zi - Gi * x[i+1],     i = n-1, n-2, ..., 1 5.2

  Blocos são: (Atil aqui é Ã) 
  - B[TB][TB]       -> TB blocos de 1 diagonal com TB elem
  - Atil[TB][3][TB] -> TB blocos de 3 diagonais de TB elem
  - G[TB][TB][TB]   -> TB blocos de TBxTB elementos
  */

  // Passo 3 >> fazendo as matrizes
  double* C = new double[TB+1];

  double** B       = new double*[TB+1];
  double** z       = new double*[TB+1];
  double** block_d = new double*[TB+1];
  double** block_w = new double*[TB+1];

  double*** Atil = new double**[TB+1];
  double*** G = new double**[TB+1];

  for(int block = 0; block <= TB; block++) {
    Atil[block] = new double*[4];
    for(int e = 0; e < 4; e++) {
      Atil[block][e] = new double[TB+1];
      for(int pos = 0; pos <= TB; pos++)
          Atil[block][e][pos] = 0.0;
    }
  }

  for(int i = 0; i <= TB; i++)
  {
    B[i]       = new double[TB+1];
    z[i]       = new double[TB+1];
    block_d[i] = new double[TB+1];
    block_w[i] = new double[TB+1];

    G[i]       = new double*[TB+1];
    for (int j = 0; j <= TB; j++)
      G[i][j] = new double[TB+1];
  }

  // Passo 3.1.1 >> Ã1 = A1
  Copy_A_to_Atil(diag, Atil[1], 1);

  // Pré: Fazer Ci
  /*
  |      [ 4 . . ] |
  | Ci = [ . 4 . ] |
  |      [ . . 4 ] |
  */
  double** C_columns = new double*[TB + 1];
  for (int i = 1; i <= TB; i++) {
    C_columns[i] = new double[TB + 1];

    for (int j = 1; j <= TB; j++)
      C_columns[i][j] = (j == i) ? diag[4][j] : 0.0; // Unit vector

    // Passo 3.1.2 >> Ã1 * G1 = C1 => G1
    LU_solve(Atil[1], G[1][i], C_columns[i]);
  }

  // Passo 3.2.1 >> Ãi = Ai - Bi * G[i-1], i = 2, 3, ..., n-1
  for(int i = 2; i < TB; i++){
    int bli = (i-1)*TB;

    for(int j = 1; j <= TB; j++) {
      for(int col = 1; col <= TB; col++)
        diag[2][bli+j] -= diag[0][bli+j] * G[i-1][j][col];
    }

    Copy_A_to_Atil(diag, Atil[i], i);

    // Passo 3.3 >> Ãi * Gi = Ci => Gi,    i = 2, 3, ..., n-1
    for (int k = 1; k <= TB; k++) {
      for (int p = 1; p <= TB; p++)
        C_columns[k][p] = (p == k) ? diag[4][bli+p] : 0.0; // Unit vector

      LU_solve(Atil[i], G[i][k], C_columns[k]);
    }
  }

  // Passo 3.2.2 >> Ãn = An - Bn * G[n-1], i = n
  for(int j = 1; j <= TB; j++)
    diag[2][DMR-TB+j] -= diag[0][DMR-TB+j] * G[TB-1][TB][j];

  Copy_A_to_Atil(diag, Atil[TB], TB);

  // Pré: fazer bloco de d1
  for(int j = 1; j <= TB; j++)
    block_d[1][j] = d[j];

  // Passo 4.1 >> Ã1 * z1 = d1 => z1
  LU_solve(Atil[1], z[1], block_d[1]);

  // Passo 4.2 >> Ãi * zi = di - Bi * z[i-1] => zi, i = 2, 3, ..., n
  for(int i = 2; i <= TB; i++){
    int bli = (i-1)*TB;
    
    // Fazendo di - Bi * z[i-1]
    for(int j = 1; j <= TB; j++)
      block_d[i][j] = d[bli+j] - diag[0][bli+j] * z[i-1][j];

    // Resolvendo Ãi * zi = di - Bi * z[i-1]
    LU_solve(Atil[i], z[i], block_d[i]);
  }

  // Passo 5.1 >> xn = zn
  for(int j = 1; j <= 25; j++)
    block_w[TB][j] = z[TB][j];

  // Passo 5.2 >> xi = zi - Gi * x[i+1],     i = n-1, n-2, ..., 1
  for(int i = TB-1; i >= 1; i--){
    for(int co = 1; co <= TB; co++){
      double sum_Gx = 0.0;
      
      for(int li = 1; li <= TB; li++)
        sum_Gx += G[i][co][li] * block_w[i+1][li];

      block_w[i][co] = z[i][co] - sum_Gx;
    }
  }

  // Passo 5.3 >> Transformar block_w to w:
  for(int i = 1; i <= TB; i++){
    int bli = (i-1)*TB;

    for(int j = 1; j <= TB; j++)
      w[bli+j] = block_w[i][j];
  }

  // Debug
  if(debug_const == 1){
    cout << fixed << setprecision(1);
    for(int i = 1; i <= TB; i++){
      for (int j = 1; j <= TB; j++)
        cout << G[2][i][j] << " ";
      cout << endl;
    }
  }

  // ============================| Print
  double Tempo_TOTAL = (double)(clock() - time_req)/CLOCKS_PER_SEC;
  double somaErro = 0.0;
  
  // Impressão dos Valores
  // Por favor, crie a pasta antes!
  ofstream file;
  string name_file = "2D_finVolMet_BLU/2D_finVolMet_BLUv2_" + to_string(DMR) + ".txt";

  file.open(name_file);

  file << "Tabela de valores para implementação com " << DMR << " subintervalos" << endl;
  file << fixed << setprecision(12);
  file << "Tempo de execução = " << Tempo_TOTAL << " seg.\n\n";
  file << "x;y;f(x,y);exact_solution;difference;Vx;Vy;VxE;VyE" << endl;

  for (int i = 1; i <= DMR; i++)
    somaErro += abs(w[i] - valor_real[i]);

  for(int i = 1; i <= M; i++){
    for(int j = 1; j <= N; j++){
      k = i + (j - 1)*M;
      file << A + (i - 0.5) * hx   << ";";
      file << Ci + (j - 0.5) * hy  << ";";
      file << w[k]                 << ";";
      file << valor_real[k]        << ";";
      file << valor_real[k] - w[k] << ";";
      file << 0.0                  << ";";
      file << 0.0                  << ";";
      file << 0.0                  << ";";
      file << 0.0                  << endl;
    }
  }

  file.close();
  cout << "(" << DMR << ") Erro médio: " << somaErro / DMR << " em " << Tempo_TOTAL << "s" << endl;


  // ============================| Cleanup
  // Deallocate memory in reverse order of allocation

  // 1. Deallocate diag (5 diagonals)
  for(int i = 0; i < 5; i++)
    delete[] diag[i];
  delete[] diag;

  // 2. Deallocate block components
  // B matrix
  for(int i = 0; i <= TB; i++)
    delete[] B[i];
  delete[] B;

  // z vectors
  for(int i = 0; i <= TB; i++)
    delete[] z[i];
  delete[] z;

  // block_d vectors
  for(int i = 0; i <= TB; i++)
    delete[] block_d[i];
  delete[] block_d;

  // block_w vectors
  for(int i = 0; i <= TB; i++)
    delete[] block_w[i];
  delete[] block_w;

  // 3. Deallocate single-dimension arrays
  delete[] C;

  // 4. Deallocate Atil 3D array (blocks -> diagonals -> elements)
  for(int block = 0; block <= TB; block++) {
    for(int diag = 0; diag < 4; diag++)
      delete[] Atil[block][diag];
    delete[] Atil[block];
  }
  delete[] Atil;

  // 5. Deallocate G 3D array (blocks -> rows -> columns)
  for(int i = 0; i <= TB; i++) {
    for(int j = 0; j <= TB; j++)
      delete[] G[i][j];
    delete[] G[i];
  }
  delete[] G;

  if(debug_const == 1)
    cout << "Cleanup complete" << endl;
}

int main()
{
  cout << "Versão 2 da Fatoração LU em Bloco" << endl;
  cout << fixed << setprecision(12);
  FinVol();

  return 0;
}

// g++ 2D_EDP_BLUv2.cpp -lm -o 2D_EDP_BLUv2 && ./2D_EDP_BLUv2