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
  
  // alpha = 0.0; // ??????
  // beta  = 1.0; // ??????

  M   = 25;
  N   = 25;

  DMR = M * N;
  
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
  hx = (B - A)/(M);
  hy = (D - C)/(N);

  hx2 = pow(hx,2);
  hy2 = pow(hy,2);

  int k;

  // Passo 1.1 >> Vetor de Ks
  Kx[0]     = K_function(A);
  Kx[DMR+1] = K_function(B);

  Ky[0]     = K_function(C);
  Ky[DMR+1] = K_function(D);

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

  for(int i = 1; i <= M; i++){
    for(int j = 1; j <= N; j++){
      k = i + (j - 1)*M;
      
      valor_real[k] = exact_solution(A + (i - 0.5) * hx, C + (j - 0.5) * hy);
    }
  }

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
      campoVel[0][k] += (i > 1) ? -(K_half(Kx[k], Kx[k - 1]) * (valor_real[k] - valor_real[k-1])) / (hx * 2) : 0.0;
      campoVel[0][k] += (i < M) ? -(K_half(Kx[k], Kx[k + 1]) * (valor_real[k+1] - valor_real[k])) / (hx * 2) : 0.0;
 
      // Passo 6.2 >> Norte + Sul
      campoVel[1][k] += (j > 1) ? -(K_half(Ky[k], Ky[k - M]) * (valor_real[k] - valor_real[k-M])) / (hy * 2) : 0.0;
      campoVel[1][k] += (j < N) ? -(K_half(Ky[k], Ky[k + 1]) * (valor_real[k+M] - valor_real[k])) / (hy * 2) : 0.0;
      
    }
  }

  // Impressão dos Valores
  // Por favor, crie a pasta antes!
  ofstream file;
  string name_file = "2D_finVolMet/true_2D_finVolMet_" + to_string(DMR) + ".txt";

  file.open(name_file);

  file << "Tabela de valores para implementação com " << DMR << " subintervalos" << endl;
  file << fixed << setprecision(12);
  file << "Tempo de execução = 0 seg.\n\n";
  file << "x;y;f(x,y);Vx;Vy" << endl;

  for(int i = 1; i <= M; i++){
    for(int j = 1; j <= N; j++){
      k = i + (j - 1)*M;
      file << A + (i - 0.5) * hx << ";";
      file << C + (j - 0.5) * hy << ";";
      file << valor_real[k]      << ";";
      file << campoVel[0][k]     << ";";
      file << campoVel[1][k]     << endl;
    }
  }

  file.close();
}

// Main
int main()
{
  cout << "Versão real" << endl;
  
  cout << fixed << setprecision(12);
  FinVol();

  return 0;
}

// g++ real_2D_EDP_CF.cpp -lm -o real_2D_EDP_CF && ./real_2D_EDP_CF