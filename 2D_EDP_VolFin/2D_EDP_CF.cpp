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

double K_function(double x) {
  return 1.0;
}

double q(double x, double y) {
  return - 8.0 * PI2 * cos(2.0 * PI * x) * cos(2.0 * PI * y);
}

// Funciona para ambos 'i' e 'j'
double K_half(double K_1, double K_2) {
  return 2 * (K_1 * K_2) / (K_1 + K_2);
}

// Função exata, se houver
double exact_solution(double x) {
  return cos(2.0 * PI * x) * cos(2.0 * PI * y);
}

// =================================================================
// Função de aproximação, com entrada de ? subintervalos.

void FinVol()
{
  // ================= | Variáveis!
  int M, N, DMR;       // Dimensões M,N da matriz e
                       // Dimensão da Matriz Resultante (DMR)

  double A, B, C, D;   // (x,y) \in [A.B] X [C,D]
  double alpha, beta;  // Soluções de contorno
  double hx, hy;       // Pulo entre x_i e x_{i+1}
  double valor_real;   // Valor real da solução
  
  int debug_const = 0; // Altere para 1 caso queira ver o tempo de
                       // execução de cada passo
 
  clock_t time_req, inter_time; // Para medir o tempo de execução!
  
  // =================================================================
  // COLOQUE AQUI OS VALORES DE INTERESSE
  
  A = C = 0.0;
  B = D = 1.0;
  
  alpha = 0.0; // ??????
  beta  = 1.0; // ??????

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
  double d[DMR+1];
  
  // RESPOSTA!
  double w[DMR + 2];
  
  // Específico
  double Kx[DMR + 2];
  double Ky[DMR + 2];

  // ============================| Começo 
  time_req = inter_time = clock();

  // Passo 1
  hx = (B - A)/(DMR + 1);
  hy = (D - C)/(DMR + 1);

  hx2 = pow(hx,2);
  hy2 = pow(hy,2);
  
  // Passo 1.1 >> Vetor de Ks
  Kx[0]     = K_function(A);
  Kx[DMR+1] = K_function(B);

  Ky[0]     = K_function(C);
  Ky[DMR+1] = K_function(D);

  for(int i = 1; i <= DMR; i++)
  {
    Kx[i] = K_function(A + (i - 0.5) * hx);
    Ky[i] = K_function(C + (i - 0.5) * hy);

    diag[2][i] = 0; // garantindo 0 para diag. principal

    // Aproveitando o loop para assinalar valores à 'd'
    // d[i] = q() >> talvez pensar numa aritmetica no futuro
  }

  if(debug_const == 1)
    cout << "Passo 1: " << (double)(clock() - inter_time)/CLOCKS_PER_SEC << endl;
  inter_time = clock();

  // Passo 2.1 >> Diagonais exteriores
  // 0: . . . 4 5 6 7 8 9 :N
  // 4: 1 2 3 4 5 6 . . . :S
  for(int i = 1; i <= DMR - M; i++)
  {
    diag[4][i]      = - K_half(Ky[i + 1]    , Ky[i]    ) / hy2; // Diag. Sul
    diag[0][M + i]  = - K_half(Ky[M + i - 1], Ky[M + i]) / hy2; // Diag. Norte

    diag[2][i]     -= diag[4][i] + diag[0][M + i];            // Diag. Central
  }

  // Passo 2.2 >> Contorno das Diagonais Exteriores
  for(int i = 1; i <= M; i++)
  {
    diag[4][DMR - i + 1] = 0; // Contorno Sul
    diag[0][i]           = 0; // Contorno Norte
  }

  if(debug_const == 1)
    cout << "Passo 2: " << (double)(clock() - inter_time)/CLOCKS_PER_SEC << endl;
  inter_time = clock();

  // Passo 3.1 >> Diagonais Internas
  /*
    OBS.: O seguinte código procura usar aritmética para acessar os valores
    das diagonais internas sem usar condicionais (if/else). O objetivo dessa
    digressão seria otimizar o programa.

    Mesmo assim, esse "pensamento otimizado" foi desenhado para um só while.
    Devido ser necessário adicionar os valores das diagonais internas à cen-
    tral, foi necessário um whileloop complementar.

    Ainda não tenho certeza se essas horas perdidas tentando encontrar uma
    aritmética ligeiramente elegante serviram de facto para otimizar o pro-
    grama, porém não gostaria de jogar fora todo o trabalho que tive por u-
    ma pedra no caminho. O código que usei como teste das diagonais, assim
    como a versão com somente um whileLoop, pode ser encontrada no arquivo
    test.cpp.

    alea iacta est.
  */

  // 1: . 2 3 . 5 6 . 8 9 :O
  int aux_geral = 0, aux_oeste = 0, aux_leste = 1;
  while(aux_oeste <= DMR)
  {
    diag[1][aux_oeste]  = - K_half(Kx[aux_oeste], Kx[aux_oeste + 1]) / hx2;

    diag[2][aux_oeste] -= diag[1][aux_oeste];

    aux_geral++;
    aux_oeste += aux_geral % 2 + 1;
  }

  // 3: 1 2 . 4 5 . 7 8 . :L
  aux_geral = 0;
  while(aux_leste <= DMR)
  {
    diag[3][aux_leste]  = - K_half(Kx[aux_leste], Kx[aux_leste + 1]) / hx2;

    diag[2][aux_leste] -= diag[3][aux_leste];

    aux_geral++;
    aux_leste += (aux_geral + 1) % 2 + 1;
  }
  
  // Passo 3.2 >> Contorno das Diagonais Internas
  int aux_contorno = 1;
  while(aux_contorno <= DMR)
  {
    diag[1][aux_contorno]   = 0;
    diag[3][aux_contorno-1] = 0;

    aux_contorno += 3;
  }
  diag[3][DMR] = 0;


  if(debug_const == 1)
    cout << "Passo 3: " << (double)(clock() - inter_time)/CLOCKS_PER_SEC << endl;
  inter_time = clock();

  // Passo 4 >> Matriz d com aritmética k (pelo livro)
  // Tudo indica que não é necessário assignment de valores 
  // fora dos forloops, como feito no exemplo 1D.

  int k; // aritmética k = i + (j - 1)M

  for(int i = 1; i <= M; i++)
  {
    for(int j = 1; j <= N; j++)
    {
      k = i + (j - 1)*M;
      d[k] = q(A + (i - 0.5) * hx, C + (j - 0.5) * hy)
    }
  }

  /*
    OBS.: Eu compreendo que depois de uma grande aritmética no
    passo 3.1, usar dois forLoops encadeados é meio broxante.
    Porém, estou cansado!
  */

  if(debug_const == 1)
    cout << "Passo 4: " << (double)(clock() - inter_time)/CLOCKS_PER_SEC << endl;
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

  // Passo 5.1 >> Início (Block I, without (Upper) bands)
  LUMatrix[0][1] = diag[0][1]; // = 0             bL
  LUMatrix[1][1] = diag[1][1]; // = 0              L
  LUMatrix[2][1] = diag[2][1];                  // C
  LUMatrix[3][1] = diag[3][1] / LUMatrix[2][1]; // U

  for (int i = 2; i <= M; i++)
  {
    LUMatrix[0][i] = diag[0][i]; // = 0
    
    LUMatrix[1][i] = diag[1][i];                                   // L
    LUMatrix[2][i] = diag[2][i] - LUMatrix[1][i]*LUMatrix[3][i-1]; // C
    LUMatrix[3][i] = diag[3][i] / LUMatrix[2][i];                  // U
  }

  // Passo 5.2 >> Intermediário (para bL, L, C e U) e Final (para bU)
  for (int i = M + 1; i <= DMR; i++)
  {
    LUMatrix[0][i] = diag[0][i];                                     // band L
    LUMatrix[4][i] = diag[4][i] / LUMatrix[2][i];                    // band U (= 0 para i > DMR - M)
    
    LUMatrix[1][i] = diag[1][i];                                     // L
    
    LUMatrix[2][i] = diag[2][i] - LUMatrix[1][i]*LUMatrix[3][i - 1]
                                - LUMatrix[0][i]*LUMatrix[4][i - M]; // C
    
    LUMatrix[3][i] = diag[3][i] / LUMatrix[2][i];                    // U (= 0 quando i = DMR)
  }








  for(int i = 2; i < N; i++)
  {
    l[i] = a[i] - (c[i] * u[i-1]);
    u[i] = b[i] / l[i];
    z[i] = (d[i] - (c[i] * z[i-1])) / l[i];
  }
  
  // Passo 7
  l[N] = a[N] - (c[N] * u[N-1]);
  z[N] = (d[N] - (c[N] * z[N-1])) / l[N];
  
  if(debug_const == 1)
    cout << "Passo 7: " << (double)(clock() - inter_time)/CLOCKS_PER_SEC << endl;

  // Passo 7
  w[0]   = alpha;
  w[N+1] = beta;
  w[N]   = z[N];

  // Passo 8
  for(int i = N - 1; i >= 1; i--)
    w[i] = z[i] - (u[i]*w[i+1]);

  double Tempo_TOTAL = (double)(clock() - time_req)/CLOCKS_PER_SEC;
  double somaErro = 0;
  
  // Impressão dos Valores
  // Por favor, crie a pasta antes!
  ofstream file;
  string name_file = "2D_finVolMet/2D_finVolMet_" + to_string(N) + ".txt";

  file.open(name_file);

  file << "Tabela de valores para implementação com " << N << " subintervalos" << endl;
  file << fixed << setprecision(12);
  file << "Tempo de execução = " << Tempo_TOTAL << " seg.\n\n";
  file << "x_i;w_i;y_i;|e_a|;" << endl;
  for (int i = 0; i <= N+1; i++)
  {
    valor_real = exact_solution(A + i*h);
    file << A + i*h << ";" << w[i] << ";" << valor_real << ";" << abs(w[i] - valor_real) << ";" << endl;
    somaErro += abs(w[i] - valor_real);
  }
  file.close();

  cout << "(" << N << ") Erro médio: " << somaErro / (N+1) << " em " << Tempo_TOTAL << "s" << endl;
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
  for(int i = 100; i < 1000; i+=100)
    FinVol(i*100);

  return 0;
}

// g++ -o EDP-CF EDP-CF.cpp && ./EDP-CF
