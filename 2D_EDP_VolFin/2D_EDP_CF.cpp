/* ================================ | Código Exercício Chemetov-Henrique | ===================================== */
//
// Tenta a implementação de Fatoração de Crout para os problemas bidimensionais do livro de Volumes Finitos.
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

// condições do lados dir e esquerdo

double K_function(double x) {
  return 1.0;
}

double q(double x, double y) {
  return 0; //4.0 * PI2 * cos(2.0 * PI * x * 0) * cos(2.0 * PI * y * 0);
}

// Funciona para ambos 'i' e 'j'
double K_half(double K_1, double K_2) {
  return 2 * (K_1 * K_2) / (K_1 + K_2);
}

// Função exata, se houver
double exact_solution(double x, double y) {
  return cos(2.0 * PI * x * 0) * cos(2.0 * PI * y * 0);
}

// =================================================================
// Função de aproximação, com entrada de ? subintervalos.

void FinVol()
{
  // ================= | Variáveis!
  int M, N, DMR;           // Dimensões M,N da matriz e
                           // Dimensão da Matriz Resultante (DMR)

  double A, B, C, D;       // (x,y) \in [A.B] X [C,D]
  double alpha, beta;      // Soluções de contorno
  double hx, hy, hx2, hy2; // Pulo entre x_i e x_{i+1}
  
  int debug_const = 0;     // Altere para 1 caso queira ver o tempo de
                           // execução de cada passo
 
  clock_t time_req, inter_time; // Para medir o tempo de execução!
  
  // =================================================================
  // COLOQUE AQUI OS VALORES DE INTERESSE
  
  A = C = 0.0;
  B = D = 1.0;
  
  M   = 25;
  N   = 25;

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

  Com indicação de ordem: 0 1 3 4    e    2
                          N O L S         C
  1: C L . S . . . . .
  2: O C L . S . . . .
  3: . O C . . S . . .
  4: N . . C L . S . .
  5: . N . O C L . S .
  6: . . N . O C . . S
  7: . . . N . . C L .
  8: . . . . N . O C L
  9: . . . . . N . O C

  1: p11 p21  .  p12  .   .   .   .   .
  2: p11 p21 p31  .  p22  .   .   .   .
  3:  .  p21 p31  .   .  p32  .   .   .
  4: p11  .   .  p12 p22  .  p13  .   .
  5:  .  p21  .  p12 p22 p32  .  p23  .
  6:  .   .  p31  .  p22 p32  .   .  p33
  7:  .   .   .  p12  .   .  p13 p23  .
  8:  .   .   .   .  p22  .  p13 p23 p33
  9:  .   .   .   .   .  p32  .  p23 p33

  diag

  0: . . . 4 5 6 7 8 9 :N
  1: . 2 3 . 5 6 . 8 9 :O
  2: 1 2 3 4 5 6 7 8 9 :C
  3: 1 2 . 4 5 . 7 8 . :L
  4: 1 2 3 4 5 6 . . . :S


  lembrando que é de baixo para cima, da esquerda para a direita
  */
  
  // Aqui é feito um primeiro array duplo para o sistema
  // penta diagonal requerido
  double** diag = new double*[5];
  for (int i = 0; i < 5; i++)
    diag[i] = new double[DMR+1];

  double x, y;
  double d[DMR+1], z[DMR+1];
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

  hx2 = pow(hx,2);
  hy2 = pow(hy,2);
  
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

      diag[2][k] = 0.0; // garantindo 0 para diag. principal

      // Aproveitando o loop para assinalar valores à 'd'
      // d[i] = q() >> talvez pensar numa aritmetica no futuro
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

      diag[0][k] = (j > 1) ? -K_half(Ky[k], Ky[k - 1]) / hy2 : 0.0;
      diag[4][k] = (j < N) ? -K_half(Ky[k], Ky[k + 1]) / hy2 : 0.0;

      // Passo 2.2 >> Diagonais internas
      // 1: . 2 3 . 5 6 . 8 9 :O
      // 3: 1 2 . 4 5 . 7 8 . :L
      diag[1][k] = (i > 1) ? -K_half(Kx[k], Kx[k - 1]) / hx2 : 0.0;
      diag[3][k] = (i < M) ? -K_half(Kx[k], Kx[k + 1]) / hx2 : 0.0;

      // Passo 2.3 >> Diagonal Central
      diag[2][k] = -(diag[0][k] + diag[1][k] + diag[3][k] + diag[4][k]);
    }
  }

  if(debug_const == 1)
    cout << "Passo 2: " << (double)(clock() - inter_time)/CLOCKS_PER_SEC << endl;
  inter_time = clock();

  // Passo 3 >> Matriz d
  for(int i = 1; i <= M; i++){
    for(int j = 1; j <= N; j++){
      k = i + (j - 1)*M;
      
      d[k]          = q(A + (i - 0.5) * hx, C + (j - 0.5) * hy);
      valor_real[k] = exact_solution(A + (i - 0.5) * hx, C + (j - 0.5) * hy);
    }
  }

  if(debug_const == 1)
    cout << "Passo 3: " << (double)(clock() - inter_time)/CLOCKS_PER_SEC << endl;
  inter_time = clock();

  // ============================| Fatoração de Crout especializada
  // Como na fatoração de Crout a matriz U tem diagonal principal igual a 1,
  // no seguinte array duplo, os índices de 0-2 compreendem à matriz L e os
  // índices de 3-4 compreendem à matriz U.
  /*
    [ 2 . . . . . . . . ]   [ _ 3 . 4 . . . . . ]
    [ 1 2 . . . . . . . ]   [ . _ 3 . 4 . . . . ]
    [ . 1 2 . . . . . . ]   [ . . _ 3 . 4 . . . ]
    [ 0 . 1 2 . . . . . ]   [ . . . _ 3 . 4 . . ]
    [ . 0 . 1 2 . . . . ] x [ . . . . _ 3 . 4 . ]
    [ . . 0 . 1 2 . . . ]   [ . . . . . _ 3 . 4 ]
    [ . . . 0 . 1 2 . . ]   [ . . . . . . _ 3 . ]
    [ . . . . 0 . 1 2 . ]   [ . . . . . . . _ 3 ]
    [ . . . . . 0 . 1 2 ]   [ . . . . . . . . _ ]

    Entende-se o numeral 1 como _ na explicação acima.

    U4: A4n / L2n -> mesmo indice
    U3: A3n / L2n -> mesmo indice
    
    L2: A2n - L0(n-3)*U4(n-3) - L1(n-1)*U3(n-1) (?)
    L1: A1n
    L0: A0n

    Ordem: L0, L1, L2, U3, U4.

    0: . . . | 4 5 6 7 8 9 :N
    1: . 2 3 | . 5 6 . 8 9 :O
    2: 1 2 3 | 4 5 6 7 8 9 :C
    3: 1 2 . | 4 5 . 7 8 . :L
    4: 1 2 3 | 4 5 6 . . . :S

    [ 2 3 .                   [       |
    [ 1 2 3   4 . . . . . ]   [   I   |       |       ]
    [ . 1 2   . 4 . . . . ]   [       |   II  |  III  ]
|             . . 4 . . . ]           |       |       ]
      [ 0 . 1 2 3 . 4 . . ]   [       |       |       ]
      [ . 0 . 1 2 3 . 4 . ] = [  IV   |   V   |   VI  ]
      [ . . 0 . 1 2 3 . 4 ]   [       |       |       ]
      [ . . . 0 . 1 2 3 . ]   [       |       |       ]
      [ . . . . 0 . 1 2 3 ]   [  VII  |  IIX  |  IX   ]
      [ . . . . . 0 . 1 2 ]   [       |       |       ]

  */

  double** LUMatrix = new double*[5];
  for (int i = 0; i < 5; i++)
    LUMatrix[i] = new double[DMR+1];

  // Passo 4.1 >> Início (Block I, without (Upper) bands)
  LUMatrix[0][1] = diag[0][1]; // = 0             bL
  LUMatrix[1][1] = diag[1][1]; // = 0              L
  LUMatrix[2][1] = diag[2][1];                  // C
  LUMatrix[3][1] = diag[3][1] / LUMatrix[2][1]; // U

  // Aproveitando loop para Lz = d
  z[1] = d[1]/LUMatrix[2][1]; // z1 = d1 / 2

  for (int i = 2; i <= M; i++)
  {
    LUMatrix[0][i] = diag[0][i]; // = 0
    
    LUMatrix[1][i] = diag[1][i];                                   // L
    LUMatrix[2][i] = diag[2][i] - LUMatrix[1][i]*LUMatrix[3][i-1]; // C
    LUMatrix[3][i] = diag[3][i] / LUMatrix[2][i];                  // U

    // zi = (di - 1.z[i-1]) / 2
    z[i] = (d[i] - LUMatrix[1][i]*z[i-1]) / LUMatrix[2][i];
  }

  // Passo 4.2 >> Intermediário (para bL, L, C e U) e Final (para bU)
  for (int i = M + 1; i <= DMR; i++)
  {
    LUMatrix[0][i] = diag[0][i];                                     // band L
    
    LUMatrix[1][i] = diag[1][i];                                     // L
    
    LUMatrix[2][i] = diag[2][i] - LUMatrix[1][i]*LUMatrix[3][i - 1]
                                - LUMatrix[0][i]*LUMatrix[4][i - M]; // C
    
    LUMatrix[4][i] = diag[4][i] / LUMatrix[2][i];                    // band U (= 0 para i > DMR - M)
    LUMatrix[3][i] = diag[3][i] / LUMatrix[2][i];                    // U (= 0 quando i = DMR)

    // zi = (di - 1.z[i-1] - 0.z[i-M]) / 2
    z[i] = (d[i] - LUMatrix[1][i]*z[i-1] - LUMatrix[0][i]*z[i-M]) / LUMatrix[2][i];
  }

  /*
    Ax = d => L(Ux) = d => Lz = d => Uw = z => w

    [ 2 . . . . . . . . ] [ a ]   [ 2a           ]   [ a =   d1               / 2 ]
    [ 1 2 . . . . . . . ] [ b ]   [ 2b + 1a      ]   [ b = ( d2 - 1.a )       / 2 ]
    [ . 1 2 . . . . . . ] [ c ]   [ 2c + 1b      ]   [ c = ( d3 - 1.b )       / 2 ]
    [ 0 . 1 2 . . . . . ] [ d ]   [ 2d + 1c + 0a ]   [ d = ( d4 - 1.c - 0.a ) / 2 ]
    [ . 0 . 1 2 . . . . ] [ e ] = [ 2e + 1d + 0b ] = [ e = ( d5 - 1.d - 0.b ) / 2 ]
    [ . . 0 . 1 2 . . . ] [ f ]   [ 2f + 1e + 0c ]   [ f = ( d6 - 1.e - 0.c ) / 2 ]
    [ . . . 0 . 1 2 . . ] [ g ]   [ 2g + 1f + 0d ]   [ g = ( d7 - 1.f - 0.d ) / 2 ]
    [ . . . . 0 . 1 2 . ] [ h ]   [ 2h + 1g + 0e ]   [ h = ( d8 - 1.g - 0.e ) / 2 ]
    [ . . . . . 0 . 1 2 ] [ i ]   [ 2i + 1h + 0f ]   [ i = ( d9 - 1.h - 0.f ) / 2 ]

    [ _ 3 . 4 . . . . . ] [ a ]   [ a + 3b + 4d ]     [ a = z1 - 3b - 4d ]
    [ . _ 3 . 4 . . . . ] [ b ]   [ b + 3c + 4e ]     [ b = z2 - 3c - 4e ]
    [ . . _ 3 . 4 . . . ] [ c ]   [ c + 3d + 4f ]     [ c = z3 - 3d - 4f ]
    [ . . . _ 3 . 4 . . ] [ d ]   [ d + 3e + 4g ]     [ d = z4 - 3e - 4g ]
    [ . . . . _ 3 . 4 . ] [ e ] = [ e + 3f + 4h ]   = [ e = z5 - 3f - 4h ]
    [ . . . . . _ 3 . 4 ] [ f ]   [ f + 3g + 4i ]     [ f = z6 - 3g - 4i ]
    [ . . . . . . _ 3 . ] [ g ]   [ g + 3h      ] M   [ g = z7 - 3h      ]
    [ . . . . . . . _ 3 ] [ h ]   [ h + 3i      ] 3   [ h = z8 - 3i      ]
    [ . . . . . . . . _ ] [ i ]   [ i           ] W   [ i = z9           ]
  */

  // Passo 5.1 >> Resolver Lz = d para obter z  
  // ...

  // Passo 5.2 >> Resolver Uw = z para obter w
  w[DMR]   = z[DMR];

  for (int i = DMR - 1; i >= DMR - M + 1; i--)
    w[i] = z[i] - LUMatrix[3][i]*w[i+1];

  for (int i = DMR - M; i >= 1; i--)
    w[i] = z[i] - LUMatrix[3][i]*w[i+1] - LUMatrix[4][i]*w[i+M];

  if(debug_const == 1)
    cout << "Passo 5: " << (double)(clock() - inter_time)/CLOCKS_PER_SEC << endl;

  // Passo 6 >> Campo de Velocidades
  double** campoVel = new double*[2];
  for (int i = 0; i < 2; i++)
    campoVel[i] = new double[DMR+1];

  for(int i = 1; i <= M; i++){
    for (int j = 1; j <= N; j++){
      k = i + (j - 1)*M;

      // Inicializando com 0
      campoVel[0][k] = 0;
      campoVel[1][k] = 0;

      // Passo 6.1 >> Leste + Oeste
      campoVel[0][k] += (i > 1) ? -(K_half(Kx[k], Kx[k - 1]) * (w[k] - w[k-1])) / (hx * 2) : 0.0;
      campoVel[0][k] += (i < M) ? -(K_half(Kx[k], Kx[k + 1]) * (w[k+1] - w[k])) / (hx * 2) : 0.0;
 
      // Passo 6.2 >> Norte + Sul
      campoVel[1][k] += (j > 1) ? -(K_half(Ky[k], Ky[k - M]) * (w[k] - w[k-M])) / (hy * 2) : 0.0;
      campoVel[1][k] += (j < N) ? -(K_half(Ky[k], Ky[k + 1]) * (w[k+M] - w[k])) / (hy * 2) : 0.0;
      
    }
  }

  // =================================================================
  // Final
  double Tempo_TOTAL = (double)(clock() - time_req)/CLOCKS_PER_SEC;
  double somaErro = 0;
  
  // Impressão dos Valores
  // Por favor, crie a pasta antes!
  ofstream file;
  string name_file = "2D_finVolMet/2D_finVolMet_" + to_string(DMR) + ".txt";

  file.open(name_file);

  file << "Tabela de valores para implementação com " << DMR << " subintervalos" << endl;
  file << fixed << setprecision(12);
  file << "Tempo de execução = " << Tempo_TOTAL << " seg.\n\n";
  file << "x;y;f(x,y);Vx;Vy" << endl;

  //int somaKappa = 0;
  for (int i = 1; i <= DMR; i++){
    //somaKappa += Ky[i];
    somaErro += abs(w[i] - valor_real[i]);
  }
  
  for(int i = 1; i <= M; i++){
    for(int j = 1; j <= N; j++){
      k = i + (j - 1)*M;
      file << A + (i - 0.5) * hx << ";";
      file << C + (j - 0.5) * hy << ";";
      file << w[k]               << ";";
      file << campoVel[0][k]     << ";";
      file << campoVel[1][k]     << endl;
    }
  }

  file.close();
  cout << "(" << DMR << ") Erro médio: " << somaErro / DMR << " em " << Tempo_TOTAL << "s" << endl;
  //cout << somaKappa / DMR << endl;

  for (int i = 0; i < 5; i++)
  {
    delete[] diag[i];
    delete[] LUMatrix[i];
  }
  delete[] diag;
  delete[] LUMatrix;
  
}

// Main
int main()
{
  cout << "Versão SEM paralelização" << endl;
  
  cout << fixed << setprecision(12);
  FinVol();

  return 0;
}

// g++ 2D_EDP_CF.cpp -lm -o 2D_EDP_CF && ./2D_EDP_CF