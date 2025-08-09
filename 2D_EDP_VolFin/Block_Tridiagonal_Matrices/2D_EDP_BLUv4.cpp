/* ================================ | Código Exercício Chemetov-Henrique | ===================================== */
//
// Tenta a implementação de Block LU factorization para matrizes em bloco tridiagonais para os problemas bidimen-
// sionais do livro de Volumes Finitos. Consultar ANALYSIS OF NUMERICAL METHODS, Isaacson e Keller.
// 
// Agora o problema é a condição de neumann que dá uma solução up to a constant. Agora vou mudar para uma condi-
// ção local dirichlet
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
// Esse problema tem condição de Dirichlet HOMOGÊNEA, dado que
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
  Atil[TB+1][TB+1][TB+1]

  .   1 2 3
  1 [ 2 3 . ]
  2 [ 1 2 3 ] = A = diag[5][DMR+1]
  3 [ . 1 2 ]

  .  1 2 3 4 5 6 7 8 9
  1: . 2 3 . 5 6 . 8 9 :O
  2: 1 2 3 4 5 6 7 8 9 :C
  3: 1 2 . 4 5 . 7 8 . :L
*/
void Copy_A_to_Atil(double **A, double **Atil, int block_index){
  int TB    = 25;
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

// Resolve sistemas lineares por decomposição LU
// Dx = f => L(Ux) = f => Lz = f => Ux = z
/*
  INPUT:
    Matriz  D[0..TB]
    Vetor   F[0..TB]
    Vetor   X[0..TB]

  Pelo que eu entendi, pelo teorema 6 (I, p.61), pos-
  so fazer decomposição LU de Ã sem pivotamento.

  O código abaixo é a fatoração LU do livro de Análi-
  se Numérica do Burden, p. 450, algoritmo 6.4.

  I: Analysis of Numerical Methods, Isaacson e Keller
*/
void LU_solve(double **D, double* X, double* F){
  int TB = 25;
  double sum;

  double* z = new double[TB+1];

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
    if (L[j][j] == 0.0) {
      throw runtime_error("Error: Singular matrix encountered in LU_solve. Division by zero.");
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

  // Passo 7.1 >> Solucionar Lz = F
  z[1] = F[1] / L[1][1];

  // Passo 7.2
  for(int i = 2; i <= TB; i++){
    sum = 0.0;
    for(int j = 1; j <= i - 1; j++)
      sum += L[i][j] * z[j];
    z[i] = (F[i] - sum)/L[i][i];
  }

  // Passo 9.1 >> Substituição regressiva
  X[TB] = z[TB];

  // Passo 9.2
  for(int i = TB - 1; i >= 1; i--){
    sum = 0.0;
    for(int j = i+1; j <= TB; j++)
      sum += U[i][j] * X[j];
    X[i] = (z[i] - sum);
  }

  // Cleanup memory
  delete[] z;
  for (int i = 0; i <= TB; i++) {
    delete[] L[i];
    delete[] U[i];
  }
  delete[] L;
  delete[] U;
}


// ====================== NEW FUNCTION TO RECONSTRUCT MATRIX A FROM DIAGONALS
double** Reconstruct_A(double** diag, int DMR, int TB) {
    double** A = new double*[DMR + 1];
    for (int i = 0; i <= DMR; i++) {
        A[i] = new double[DMR + 1];
        for (int j = 0; j <= DMR; j++) A[i][j] = 0.0;
    }

    for (int i = 1; i <= DMR; i++) {
        // Lower diagonal (y-direction)
        if (i > TB) A[i][i - TB] = diag[0][i];
        // Lower diagonal (x-direction)
        if (i > 1) A[i][i - 1] = diag[1][i];
        // Main diagonal
        A[i][i] = diag[2][i];
        // Upper diagonal (x-direction)
        if (i < DMR) A[i][i + 1] = diag[3][i];
        // Upper diagonal (y-direction)
        if (i <= DMR - TB) A[i][i + TB] = diag[4][i];
    }
    return A;
}

void FinVol(int debug_const)
{
  // ================= | Variáveis!
  int M, N; // Dimensões M,N da matriz e
  int DMR;  // Dimensão da Matriz Resultante (DMR)
  int TB;   // Tamanho do bloco

  double Al, Be, Ga, De;   // (x,y) \in [Al.Be] X [Ga,De]
  double hx, hy, hx2, hy2; // Pulo entre x_i e x_{i+1}
  
  // int debug_const = 0;  // Altere para 1 caso queira ver o tempo de
                           // execução de cada passo
 
  clock_t time_req, inter_time; // Para medir o tempo de execução!
  
  // =================================================================
  // COLOQUE AQUI OS VALORES DE INTERESSE
  
  Al = Ga = 0.0;
  Be = De = 1.0;
  
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
  1: 2 3 . 4 . . . . .   p11       d1   k = 1 + (1 - 1).3 = 1 + 0 = 1
  2: 1 2 3 . 4 . . . .   p21       d2   k = 2 + (1 - 1).3 = 2 + 0 = 2
  3: . 1 2 . . 4 . . .   p31       d3   k = 3 + (1 - 1).3 = 3 + 0 = 3
  4: 0 . . 2 3 . 4 . .   p12       d4   k = 1 + (2 - 1).3 = 1 + 3 = 4
  5: . 0 . 1 2 3 . 4 .   p22   =   d5   k = 2 + (2 - 1).3 = 2 + 3 = 5
  6: . . 0 . 1 2 . . 4   p32       d6   k = 3 + (2 - 1).3 = 3 + 3 = 6
  7: . . . 0 . . 2 3 .   p13       d7   k = 1 + (3 - 1).3 = 1 + 6 = 7
  8: . . . . 0 . 1 2 3   p23       d8   k = 2 + (3 - 1).3 = 2 + 6 = 8
  9: . . . . . 0 . 1 2   p33       d9   k = 3 + (3 - 1).3 = 3 + 6 = 9
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
  hx = (Be - Al)/M;
  hy = (De - Ga)/N;

  hx2 = hx * hx;
  hy2 = hy * hy;
  
  // Passo 1.1 >> Vetor de Ks
  Kx[0]     = K_function(Al);
  Kx[DMR+1] = K_function(Be);

  Ky[0]     = K_function(Ga);
  Ky[DMR+1] = K_function(De);

  for(int i = 1; i <= M; i++){
    for(int j = 1; j <= N; j++){
      k = i + (j - 1)*M;

      Kx[k] = 1.0;
      Ky[k] = 1.0;

      // Passo 1.2 >> Vetor 'd' e das soluções exatas 
      d[k]          = q(Al + (i - 0.5) * hx, Ga + (j - 0.5) * hy);
      valor_real[k] = exact_solution(Al + (i - 0.5) * hx, Ga + (j - 0.5) * hy);
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
      diag[0][k] = (j > 1) ? -1.0 / hy2 : 0.0;
      diag[4][k] = (j < N) ? -1.0 / hy2 : 0.0; 

      // Passo 2.2 >> Diagonais internas (x-direction)
      diag[1][k] = (i > 1) ? -1.0 / hx2 : 0.0;
      diag[3][k] = (i < M) ? -1.0 / hx2 : 0.0;

      // Passo 2.3 >> Diagonal Central
      diag[2][k] = - diag[0][k] - diag[1][k] - diag[3][k] - diag[4][k];
    }
  }

  // Passo 2.4 >> Introduzindo solução de Dirichlet local em (1,1)

  int k_ref = 1;

  diag[2][k_ref] = 1.0;                  // Diagonal principal igual a 1
  diag[4][k_ref] = diag[3][k_ref] = 0.0; // Resto zerado
  
  d[k_ref] = valor_real[k_ref];

  if(debug_const == 1){
    if(fabs(diag[2][1]) <= fabs(diag[3][1]))
      cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    for(int k = 1; k < DMR; k++){
      if(fabs(diag[2][k]) < fabs(diag[1][k]) + fabs(diag[3][k]))
        cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    }
    if(fabs(diag[2][DMR]) <= fabs(diag[1][DMR]))
      cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
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
    [ B2 Ã2    ] [    I2 G2 ] = LU
    [    B3 Ã3 ] [       I3 ]

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

  Blocos são:
  - B   [TB][TB][TB] -> TB blocos de TBxTB elementos
  - Atil[TB][TB][TB] -> ""
  - G   [TB][TB][TB] -> ""

  - C   [TB][TB]     -> TBxTB elementos
  - Z   [TB][TB]     -> ""

  Para debug, vou criar matrizes Matrix, L e U de veri-
  ficação da Fatoração LU em bloco.
  */

  // DEBUG
  double**  Matrix = new double*[DMR+1];
  double**  L      = new double*[DMR+1];
  double**  U      = new double*[DMR+1];
  for(int i = 0; i <= DMR; i++){
    Matrix[i] = new double[DMR+1];
    L[i]      = new double[DMR+1];
    U[i]      = new double[DMR+1];
    for(int j = 0; j <= DMR; j++){
      Matrix[i][j] = 0.0;
      L[i][j]      = 0.0;
      U[i][j]      = 0.0;
    }
  }

  // Passo 3 >> fazendo as matrizes
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

    // Alocando linhas
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

  // Passo 3.1.1 >> Ã1 = A1
  Copy_A_to_Atil(diag, Atil[1], 1);

  // Passo 3.1.2 >> Ã1 * G1 = C1 => G1
  for(int i = 1; i <= TB; i++){
    // Pré: fazer C1 (vou ter q fazer isso sempre)
    C[i][i] = diag[4][i];

    LU_solve(Atil[1], G[1][i], C[i]);
  }

  // Passo 3.2.1 >> Ãi = Ai - Bi * G[i-1], i = 2, 3, ..., n-1
  for(int b = 2; b < TB; b++){
    int bli = (b-1)*TB;

    // Primeiro copio os valores de diag para Atil no bloco
    Copy_A_to_Atil(diag, Atil[b], b);

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

      LU_solve(Atil[b], G[b][i], C[i]);
    }
  }

  // Passo 3.2.2 >> Ãn = An - Bn * G[n-1], i = n
  int bli = (TB-1)*TB;

  Copy_A_to_Atil(diag, Atil[TB], TB);

  for(int i = 1; i <= TB; i++){
    for(int j = 1; j <= TB; j++){
      Atil[TB][i][j] -= diag[0][bli+i] * G[TB-1][i][j];
    }
  }

  if(debug_const == 1)
    cout << "Passo 3: " << (double)(clock() - inter_time)/CLOCKS_PER_SEC << endl;
  inter_time = clock();

  // Pré: fazer bloco de d
  for(int b = 1; b <= TB; b++){
    int bli = (b-1)*TB;
    for(int j = 1; j <= TB; j++)
      b_d[b][j] = d[bli + j];
  }

  // Passo 4.1 >> Ã1 * z1 = d1 => z1
  LU_solve(Atil[1], Z[1], b_d[1]);

  // Passo 4.2 >> Ãi * zi = di - Bi * z[i-1] => zi, i = 2, 3, ..., n
  for(int b = 2; b <= TB; b++){
    if(debug_const == 1){
      cout << fixed << setprecision(12);
      cout << "b" << b << ": ";
    }

    int bli = (b-1)*TB;
    // Fazendo di - Bi * z[i-1]
    for(int j = 1; j <= TB; j++){
      b_i[j] = b_d[b][j] - diag[0][bli+j] * Z[b-1][j];
      if(debug_const == 1) cout << b_i[j] << " ";
    }
    if(debug_const == 1) cout << endl;

    // Resolvendo Ãi * zi = bi = di - Bi * z[i-1]
    LU_solve(Atil[b], Z[b], b_i);
  }

  if(debug_const == 1){
    cout << fixed << setprecision(1);
    int bli = (TB-1)*TB;

    for(int i = 1; i <= TB; i++){
      cout << i << " ";
      cout << diag[0][i] << " ";
      cout << diag[1][i] << " ";
      cout << diag[2][i] << " ";
      cout << diag[3][i] << " ";
      cout << diag[4][i] << endl;
    }
  }
  

  if(debug_const == 1){
    cout << fixed << setprecision(4);
    
    for(int b = 1; b <= TB; b++){
      cout << "Z" << b <<": ";
      for (int j = 1; j <= TB; j++)
        cout << Z[b][j] << " ";
      cout << endl;

    }
  }

  if(debug_const == 1)
    cout << "Passo 4: " << (double)(clock() - inter_time)/CLOCKS_PER_SEC << endl;
  inter_time = clock();

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

  if(debug_const == 1){
    cout << fixed << setprecision(1);
    for(int b = 1; b <= TB; b++){
      cout << "W" << b <<": ";
      for (int j = 1; j <= TB; j++)
        cout << b_w[b][j] << " ";
      cout << endl;
    }
  }

  if(debug_const == 1){
    cout << fixed << setprecision(1);
    for(int b = 1; b <= TB; b++){
      cout << "d" << b <<": ";
      for (int j = 1; j <= TB; j++)
        cout << b_d[b][j] << " ";
      cout << endl;
    }
  }

  if(debug_const == 1){
    cout << fixed << setprecision(10);
    for(int i = 1; i <= TB; i++){
      for(int j = 1; j <= TB; j++)
        cout << Atil[25][i][j] << " ";
      cout << endl;
    }
  }

  // Passo 5.3 >> Transformar block_w to w:
  for(int i = 1; i <= TB; i++){
    int bli = (i-1)*TB;

    for(int j = 1; j <= TB; j++)
      w[bli+j] = b_w[i][j];
  }

  // ============================| Print
  double Tempo_TOTAL = (double)(clock() - time_req)/CLOCKS_PER_SEC;
  double somaErro = 0.0;
  
  // Impressão dos Valores
  // Por favor, crie a pasta antes!
  ofstream file;
  string name_file = "2D_finVolMet_BLU/2D_finVolMet_BLUv4_" + to_string(DMR) + ".txt";

  file.open(name_file);

  file << "Tabela de valores para implementação com " << DMR << " subintervalos" << endl;
  file << fixed << setprecision(12);
  file << "Tempo de execução = " << Tempo_TOTAL << " seg.\n\n";
  file << "x;y;f(x,y);exact_solution;difference;Vx;Vy;VxE;VyE" << endl;

  for (int i = 1; i <= DMR; i++)
    somaErro += abs(w[i] - valor_real[i]);

  file << fixed << setprecision(12);
  for(int i = 1; i <= M; i++){
    for(int j = 1; j <= N; j++){
      k = i + (j - 1)*M;                                      // Erro aqui?                   <=======================
      file << Al + (i - 0.5) * hx  << ";";
      file << Ga + (j - 0.5) * hy  << ";";
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

  // ============================| DEBUG
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

  [ Ã1       ] [ I1 G1    ]
  [ B2 Ã2    ] [    I2 G2 ] = LU
  [    B3 Ã3 ] [       I3 ]
  */

  // Passo D.1 >> L
  // Ã1
  for(int i = 1; i <= TB; i++){
    for(int j = 1; j <= TB; j++)
      L[i][j] = Atil[1][i][j];
  }

  // Ãb
  for(int b = 2; b <= TB; b++){
    int bli = (b-1)*TB;

    for(int i = 1; i <= TB; i++){
      for(int j = 1; j <= TB; j++)
        L[bli + i][bli + j] = Atil[b][i][j];
    }
  }

  // Populate B_i blocks in L (lower block diagonal)
  for (int b = 2; b <= TB; b++) {
    int row_start = (b - 1) * TB + 1;
    int col_start = (b - 2) * TB + 1;
    for (int i = 1; i <= TB; i++) {
      L[row_start + i - 1][col_start + i - 1] = diag[0][row_start + i - 1];
    }
  }

  // Passo D.2 >> U
  // Ii
  for(int i = 1; i <= DMR; i++)
    U[i][i] = 1.0;

  // Gi
  for (int b = 1; b < TB; b++) {
    int bli1 = b * TB;      // Start index of the next block
    int bli2 = (b - 1) * TB; // Start index of the current block
    for (int i = 1; i <= TB; i++) {
      for (int j = 1; j <= TB; j++) {
        U[bli2 + i][bli1 + j] = G[b][i][j]; // G placed in upper off-diagonal
      }
    }
  }

  // Passo D.3 >> Compute M = L * U with proper matrix multiplication
  for (int i = 1; i <= DMR; i++) {
    for (int j = 1; j <= DMR; j++) {
      Matrix[i][j] = 0.0;
      for (int k = 1; k <= DMR; k++) {
        Matrix[i][j] += L[i][k] * U[k][j];
      }
    }
  }

  // ====================== VALIDATION: Compare M with original A ======================
  double** A_original = Reconstruct_A(diag, DMR, TB);
  double max_error = 0.0;
  for (int i = 1; i <= DMR; i++) {
      for (int j = 1; j <= DMR; j++) {
          double diff = fabs(Matrix[i][j] - A_original[i][j]);
          if (diff > max_error) max_error = diff;
      }
  }
  cout << "VALIDATION: Max error between M (L*U) and A: " << scientific << max_error << endl;

  // Impressão dos Valores
  // Por favor, crie a pasta antes!
  ofstream debug_file;
  string name_debug_file = "2D_finVolMet_BLU/debugMatrix.txt";

  debug_file.open(name_debug_file);

  debug_file << "x;y;M;L;U" << endl;

  debug_file << fixed << setprecision(12);
  for(int i = 1; i <= DMR; i++){
    for(int j = 1; j <= DMR; j++){
      debug_file << i            << ";";
      debug_file << j            << ";";
      debug_file << Matrix[i][j] << ";";
      debug_file << L[i][j]      << ";";
      debug_file << U[i][j]      << endl;
    }
  }

  // ============================| Cleanup

  // Cleanup A_original
  for (int i = 0; i <= DMR; i++) delete[] A_original[i];
  delete[] A_original;

  // Cleanup diag array
  for(int i = 0; i < 5; i++)
      delete[] diag[i];
  delete[] diag;

  // Cleanup Atil, G, and B arrays
  for(int i = 0; i <= TB; i++) {

    // Cleanup Atil[i], G[i], B[i]
    if (Atil[i] != nullptr) {
      for(int j = 0; j <= TB; j++)
        delete[] Atil[i][j];
      delete[] Atil[i];
    }
    if (G[i] != nullptr) {
      for(int j = 0; j <= TB; j++)
        delete[] G[i][j];
      delete[] G[i];
    }
    if (B[i] != nullptr) {
      for(int j = 0; j <= TB; j++)
        delete[] B[i][j];
      delete[] B[i];
    }

    // Cleanup C[i], Z[i], b_d[i], b_w[i]
    delete[] C[i];
    delete[] Z[i];
    delete[] b_d[i];
    delete[] b_w[i];
  }

  // Delete the top-level arrays
  delete[] Atil;
  delete[] G;
  delete[] B;
  delete[] C;
  delete[] Z;
  delete[] b_d;
  delete[] b_w;
  delete[] b_i;

  for(int i = 0; i <= DMR; i++){
    delete[] Matrix[i];
    delete[] L[i];
    delete[] U[i];
  }
  delete[] Matrix;
  delete[] L;  
  delete[] U;

  if(debug_const == 1)
    cout << "Cleanup complete" << endl;
}

int main()
{
  cout << "Versão 4 da Fatoração LU em Bloco" << endl;
  cout << fixed << setprecision(12);
  FinVol(0);

  return 0;
}

// g++ 2D_EDP_BLUv4.cpp -lm -o 2D_EDP_BLUv4 && ./2D_EDP_BLUv4