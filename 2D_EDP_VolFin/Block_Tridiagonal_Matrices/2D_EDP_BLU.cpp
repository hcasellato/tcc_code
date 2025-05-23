/* ================================ | Código Exercício Chemetov-Henrique | ===================================== */
//
// Tenta a implementação de Block LU factorization para matrizes em bloco tridiagonais para os problemas bidimen-
// sionais do livro de Volumes Finitos. Consultar ANALYSIS OF NUMERICAL METHODS, Isaacson e Keller.
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

int can_i_delete_now = 0;

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

// Decompõe e resolve matrizes tridiagonais por decomposição LU
// Dx = f => L(Ux) = f => Lz = f => Ux = y
/*
  INPUT:
    Matriz  (D) das diagonais,
    Início (in) da matriz,
    Vetor   (F),
    Vetor   (x).
*/

// Passo 1 >> Alocação de vetores doubles
double* l = new double[25 + 1];  // l[1..T]
double* u = new double[25 + 1];  // u[1..T-1]
double* z = new double[25 + 1];  // z[1..T]

void LU_solve(double **D, double* F, double* x, int in) {
  /*

  [ D21 D31  .   .   .  ] = [ l_1  .   .   .   .  ] [  1 u_1  .   .   .  ]
  [ D12 D22 D32  .   .  ] = [ D12 l_2  .   .   .  ] [  .  1  u_2  .   .  ]
  [  .  D13 D23 D33  .  ] = [  .  D13 l_3  .   .  ] [  .   .  1  u_3  .  ]
  [  .   .  D14 D24 D34 ] = [  .   .  D14 l_4  .  ] [  .   .   .  1  u_4 ]
  [  .   .   .  D15 D25 ] = [  .   .   .  D15 l_5 ] [  .   .   .   .  1  ]
  
  */
  int T = 25;
  int fi = in + T - 1;
  int gi;
  
  l[0] = 0;
  u[0] = 0;
  z[0] = 0; 

  // Passo 2.1 >> Passos para decomposição LU
  l[1] = D[2][in];
  u[1] = D[3][in] / l[1];
  z[1] = F[in] / l[1];
  
  // Passo 2.2
  for(int i = 2; i < T; i++)
  {
    gi = in + i - 1;

    l[i] = D[2][gi] - (D[1][gi] * u[i-1]);
    u[i] = D[3][gi] / l[i];
    z[i] = (F[gi] - (D[1][gi] * z[i-1])) / l[i];
  }

  // Passo 2.3
  gi = in + T - 1;
  l[T] = D[2][gi] - (D[1][gi] * u[T-1]);
  z[T] = (F[gi] - (D[1][gi] * z[T-1])) / l[T];
  
  // Passo 3.1 >> Passos para 'backward substitution'
  x[T] = z[T];

  // Passo 3.2
  for(int i = T - 1; i >= 1; i--)
    x[i] = z[i] - u[i] * x[1 + i];

}

void FinVol()
{
  // ================= | Variáveis!
  int M, N; // Dimensões M,N da matriz e
  int DMR;  // Dimensão da Matriz Resultante (DMR)
  int TB;   // Tamanho do bloco

  double A, B, C, D;       // (x,y) \in [A.B] X [C,D]
  double alpha, beta;      // Soluções de contorno
  double hx, hy, hx2, hy2; // Pulo entre x_i e x_{i+1}
  
  int debug_const = 1;     // Altere para 1 caso queira ver o tempo de
                           // execução de cada passo
 
  clock_t time_req, inter_time; // Para medir o tempo de execução!
  
  // =================================================================
  // COLOQUE AQUI OS VALORES DE INTERESSE
  
  A = C = 0.0;
  B = D = 1.0;
  
  M   = 25;
  N   = 25;

  DMR = M * N;
  TB  = 25; // (int)sqrt((double)DMR);
  // if (TB * TB < DMR) TB++; I really tried ;~;

  int bli = 0; // Block Local Index

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
  
  // Aqui é feito um primeiro array duplo para o sistema
  // penta diagonal requerido
  double** diag = new double*[5];
  for(int i = 0; i < 5; i++)
    diag[i] = new double[DMR+1];

  double x, y;
  double d[DMR+1];
  double valor_real[DMR+1];   // Vetor de valores reais da solução
  
  // RESPOSTA!
  double w[DMR + 2];
  
  // Específico
  double Kx[DMR + 2];
  double Ky[DMR + 2];

  // ============================| Começo 
  time_req = inter_time = clock();

  // Passo 1
  hx = (B - A)/M;
  hy = (D - C)/N;

  hx2 = hx * hx;
  hy2 = hy * hy;
  
  // Passo 1.1 >> Vetor de Ks
  Kx[0]     = K_function(A);
  Kx[DMR+1] = K_function(B);

  Ky[0]     = K_function(C);
  Ky[DMR+1] = K_function(D);

  int k; // aritmética k = i + (j - 1)M

  for(int i = 1; i <= M; i++)
  {
    for(int j = 1; j <= N; j++)
    {
      k = i + (j - 1)*M;

      Kx[k] = K_function(A + (i - 0.5) * hx);
      Ky[k] = K_function(C + (j - 0.5) * hy);

      // Passo 1.2 >> Vetor 'd' e das soluções exatas 
      d[k]          = q(A + (i - 0.5) * hx, C + (j - 0.5) * hy);
      valor_real[k] = exact_solution(A + (i - 0.5) * hx, C + (j - 0.5) * hy);

    }
  }

  if(debug_const == 1)
    cout << "Passo 1: " << (double)(clock() - inter_time)/CLOCKS_PER_SEC << endl;
  inter_time = clock();

  // Passo 2 >> Montagem da Matriz
  for(int i = 1; i <= M; i++){
    for (int j = 1; j <= N; j++){
      k = i + (j - 1)*M;

      // Passo 2.1 >> Diagonais exteriores
      // 0: . . . 4 5 6 7 8 9 :N
      // 4: 1 2 3 4 5 6 . . . :S
      diag[0][k] = (j > 1) ? -(1.0) / hy2 : 0.0;
      diag[4][k] = (j < N) ? -(1.0) / hy2 : 0.0;

      // Passo 2.2 >> Diagonais internas
      // 1: . 2 3 . 5 6 . 8 9 :O
      // 3: 1 2 . 4 5 . 7 8 . :L
      diag[1][k] = (i > 1) ? -(1.0) / hx2 : 0.0;
      diag[3][k] = (i < M) ? -(1.0) / hx2 : 0.0;

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
    
    Ã1 = A1, G1 = A1` * C1 (.` é a inversa)        3.1
    
    Ãi = Ai - B[i] * G[i-1], i = 2, 3, ..., n      3.2
    Gi = Ai` * Ci,           i = 2, 3, ..., n-1    3.3
  
    obs1.: Ai * Gi = Ai * Ai` * Ci => Ai * Gi = Ci

    então,

    y1 = Ã1` * f1                                   4.1
    yi = Ãi` * (fi - Bi * y[i-1]), i = 2, 3, ..., n 4.2
    
    obs2.: Ã1 * y1 = Ã1 * Ã1` * f1 => Ã1 * y1 = f1
    
    xn = yn                                             5.1
    xi = yi - Gi * x[i+1],         i = n-1, n-2, ..., 1 5.2

  */

  double** LUMatrix = new double*[5];
  for(int i = 0; i < 5; i++)
    LUMatrix[i] = new double[DMR+1];

  double** gamma = new double*[TB+1]; // Vetor de gamma com TB blocos
  double** z     = new double*[TB+1]; // Vetor intermediario Lz = f
  
  for(int i = 0; i <= TB; i++)
  {
    gamma[i] = new double[TB+1];
    z[i]     = new double[TB+1];
  }

  /*
  Explicação para LUMatrix:

    LUMatrix[4][1] = Upper Band         => Sul

    LUMatrix[3][1] = Diagonal Superior  => Leste
    LUMatrix[2][1] = Diagonal Principal => Centro
    LUMatrix[1][1] = Diagonal Inferior  => Oeste
    
    LUMatrix[0][1] = Lower Band         => Norte

  Relembrando:
  [ A1 C1    ]   [ Ã1       ] [ I1 G1    ]
  [ B1 A2 C2 ] = [ B1 Ã2    ] [    I2 G2 ] = LU
  [    B2 A3 ]   [    B2 Ã3 ] [       I3 ]

|      [ 2 3 . ]       [ 4 . . ]       [ 0 . . ]
| Ai = [ 1 2 3 ], Ci = [ . 4 . ], Bi = [ . 0 . ]
|      [ . 1 2 ]       [ . . 4 ]       [ . . 0 ]

  */

  // eu acho que não precisa mais usar k, pq estou 
  // fazendo a solução de LU, certo?
  
  // Passo 3.1 >> Início (Cópia de tudo)
  for(int k = 1; k <= DMR; k++){
    // Ã1 = A1
    LUMatrix[3][k] = diag[3][k];
    LUMatrix[2][k] = diag[2][k];
    LUMatrix[1][k] = diag[1][k];
  }

  // A1 * G1 = C1 => G1
  LU_solve(diag, diag[4], gamma[1], 1);

  int inicio, final;
  for(int r = 2; r < TB; r++)
  {
    inicio = (r-1)*TB + 1;
    final  = r*TB;

    // Passo 3.2.1
    for(int k = inicio; k <= final; k++)
    {
      bli = k - inicio + 1;
      // B[i]
      LUMatrix[0][k] = diag[0][k];

      // Ãi = A[i] - B[i] * G[i-1], i = 2, 3, ..., n-1
      LUMatrix[2][k] = diag[2][k] - LUMatrix[0][k] * gamma[r-1][bli];

      // É só a diagonal central, pois B e Gamma só tem diagonal central
    }

    // Passo 3.3
    // Ai * Gi = Ci,           i = 2, 3, ..., n-1 
    LU_solve(diag, diag[4], gamma[r], inicio);
  }


  // Passo 3.2.2
  for(int k = DMR - TB + 1; k <= DMR; k++)
  {
    bli = k - (DMR - TB + 1) + 1;
    // B[i]
    LUMatrix[0][k] = diag[0][k];

    // Ãi = A[i] - B[i] * G[i-1], i = n 
    LUMatrix[2][k] = diag[2][k] - LUMatrix[0][k] * gamma[TB][bli];
  }

  int gi = 1;
  for(int r = 1; r < TB; r++){
    for(int k = 1; k <= TB; k++){
      LUMatrix[4][gi] = gamma[r][k];
      gi++;
    }
  }

  if(debug_const == 1)
    cout << "Passo 3: " << (double)(clock() - inter_time)/CLOCKS_PER_SEC << endl;
  inter_time = clock();

  // Passo 4.1
  // Ã1 * y1 = f1 => y1
  LU_solve(LUMatrix, d, z[1], 1);

  // intermediate vector for LU_Solve
  // double* i_vector = new double[TB + 1];

  // Passo 4.2
  // Ã1 * yi = fi - Bi * y[i-1], i = 2, 3, ..., n
  for(int r = 2; r <= TB; r++)
  {
    inicio = (r-1)*TB + 1;
    final  = r*TB;
    for(int k = inicio; k <= final; k++)
      d[k] -= LUMatrix[0][k] * z[r-1][k - inicio + 1];

    LU_solve(LUMatrix, d, z[r], inicio);
  }
  if(debug_const == 1)
    cout << "Passo 4: " << (double)(clock() - inter_time)/CLOCKS_PER_SEC << endl;
  inter_time = clock();

  // Passo 5.1
  // xn = yn
  for(int k = DMR; k >= DMR - TB + 1; k--){
    bli  = k - (DMR - TB + 1) + 1;
    w[k] = z[TB][bli];
  }

  // Passo 5.2
  // xi = yi - Gi * x[i+1], i = n-1, n-2, ..., 1
  for(int r = TB - 1; r >= 1; r--)
  {
    inicio = r*TB;
    final  = (r-1)*TB + 1;
    for(int k = inicio; k >= final; k--)
      w[k] = z[r][k - final + 1] - gamma[r][k - final + 1] * w[k + TB];
  }


  if(debug_const == 1)
    cout << "Passo 5: " << (double)(clock() - inter_time)/CLOCKS_PER_SEC << endl;
  inter_time = clock();

  // Passo 6 >> Campo de Velocidades
  double** campoVel = new double*[4];
  for(int i = 0; i < 4; i++)
    campoVel[i] = new double[DMR+1];

  for(int i = 1; i <= M; i++){
    for(int j = 1; j <= N; j++){
      k = i + (j - 1)*M;

      // Inicializando com 0
      campoVel[0][k] = 0;
      campoVel[1][k] = 0;

      // Passo 6.1 >> Leste + Oeste
      campoVel[0][k] += (i > 1) ? -((1.0) * (w[k] - w[k-1])) / (hx * 2) : 0.0;
      campoVel[0][k] += (i < M) ? -((1.0) * (w[k+1] - w[k])) / (hx * 2) : 0.0;
 
      campoVel[2][k] += (i > 1) ? -((1.0) * (valor_real[k] - valor_real[k-1])) / (hx * 2) : 0.0;
      campoVel[2][k] += (i < M) ? -((1.0) * (valor_real[k+1] - valor_real[k])) / (hx * 2) : 0.0;

      // Passo 6.2 >> Norte + Sul
      campoVel[1][k] += (j > 1) ? -((1.0) * (w[k] - w[k-M])) / (hy * 2) : 0.0;
      campoVel[1][k] += (j < N) ? -((1.0) * (w[k+M] - w[k])) / (hy * 2) : 0.0;

      campoVel[3][k] += (j > 1) ? -((1.0) * (valor_real[k] - valor_real[k-M])) / (hy * 2) : 0.0;
      campoVel[3][k] += (j < N) ? -((1.0) * (valor_real[k+M] - valor_real[k])) / (hy * 2) : 0.0;    }
  }

  if(debug_const == 1)
    cout << "Passo 6: " << (double)(clock() - inter_time)/CLOCKS_PER_SEC << endl;
  inter_time = clock();

  // =================================================================
  // Final
  double Tempo_TOTAL = (double)(clock() - time_req)/CLOCKS_PER_SEC;
  double somaErro = 0;
  
  // Impressão dos Valores
  // Por favor, crie a pasta antes!
  ofstream file;
  string name_file = "2D_finVolMet_BLU/2D_finVolMet_BLU_" + to_string(DMR) + ".txt";

  file.open(name_file);

  file << "Tabela de valores para implementação com " << DMR << " subintervalos" << endl;
  file << fixed << setprecision(12);
  file << "Tempo de execução = " << Tempo_TOTAL << " seg.\n\n";
  file << "x;y;f(x,y);exact_solution;difference;Vx;Vy;VxE;VyE" << endl;

  //int somaKappa = 0;
  // cout << fixed << setprecision(1);
  // cout << "N S O L C" << endl;
  for (int i = 1; i <= DMR; i++){
    cout << i << " ";
    cout << LUMatrix[0][i] << " ";
    cout << LUMatrix[1][i] << " ";
    cout << LUMatrix[2][i] << " ";
    cout << LUMatrix[3][i] << " ";
    cout << LUMatrix[4][i] << " ";
    cout << w[i]           << endl;
    somaErro += abs(w[i] - valor_real[i]);
  }
  cout << "====================\nr k gi" << endl;

  gi = 1;
  for(int r = 1; r <= TB; r++){
    for(int k = 1; k <= TB; k++){
      cout << r << " " << k << " " << gi << " ";
      cout << gamma[r][k] << " " << z[r][k] << endl;
      gi++;
    }
  }

  for(int i = 1; i <= M; i++){
    for(int j = 1; j <= N; j++){
      k = i + (j - 1)*M;
      file << A + (i - 0.5) * hx   << ";";
      file << C + (j - 0.5) * hy   << ";";
      file << w[k]                 << ";";
      file << valor_real[k]        << ";";
      file << valor_real[k] - w[k] << ";";
      file << campoVel[0][k]       << ";";
      file << campoVel[1][k]       << ";";
      file << campoVel[2][k]       << ";";
      file << campoVel[3][k]       << endl;
    }
  }

  file.close();
  cout << "(" << DMR << ") Erro médio: " << somaErro / DMR << " em " << Tempo_TOTAL << "s" << endl;
  //cout << somaKappa / DMR << endl;

  for(int i = 0; i < 5; i++)
  {
    delete[] diag[i];
    delete[] LUMatrix[i];
  }
  delete[] diag;
  delete[] LUMatrix;

  for (int i = 1; i <= TB; i++) {
    delete[] gamma[i];
    delete[] z[i];
  }
  delete[] gamma;
  delete[] z;

  for (int i = 0; i < 4; i++)
    delete[] campoVel[i];
  delete[] campoVel;
}

// Main
int main()
{
  cout << "Versão SEM paralelização" << endl;
  
  cout << fixed << setprecision(12);
  FinVol();

  delete[] l;
  delete[] u;
  delete[] z;


  return 0;
}

// g++ 2D_EDP_BLU.cpp -lm -o 2D_EDP_BLU && ./2D_EDP_BLU