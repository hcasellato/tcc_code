/* ================================ | Código Exercício Chemetov-Henrique | ===================================== */
//
// Tenta a implementação de Fatoração de Crout para os problemas unidimensionais do livro de Volumes Finitos.
//
// Problemas da forma: 
//
// - d/dx (K dp/dx) = q
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

// =================================================================
// PROBLEMA:
//
// - d/dx (K dp/dx) = -25cos(25x) em OMEGA = [0,1]
//                p = x           na borda de OMEGA
//
// Onde K(x) = 2 + sen(25x)
// Sol. Exata p(x) = x
// =================================================================

double K_function(double x) {
  return 2.0 + sin(25.0 * x);
}

double q(double x) {
  return - 25.0 * cos(25.0 * x);
}

double K_half(double K_1, double K_2) {
  // Pode ser pensada em um array de Ks, mas
  // depende da futura implementação
  
  return 2 * (K_1 * K_2) / (K_1 + K_2);
}

// Função exata, se houver
double exact_solution(double x) {

    return x;
}

// =================================================================
// Função de aproximação, com entrada de N subintervalos.
// Ainda há de ser generalizada, porém é útil, por ora, assim mesmo

void FinVol(int N)
{
  // ================= | Variáveis!
  double A, B;         // Limites dos Intervalos e.g. x \in [A.B]
  double alpha, beta;  // Soluções de contorno
  double h;
  double valor_real;
  
  int debug_const = 0; // altere para 1 caso queira ver o tempo de
                       // exec de cada passo

  //int N; // Número de Subintervalos
 
  clock_t time_req, inter_time; // Para medir o tempo de execução!
  
  // =================================================================
  // COLOQUE AQUI OS VALORES DE INTERESSE
  
  // N = 19;    // Número de Subintervalos
  A = 0.0;
  B = 1.0;
  
  alpha = 0.0; // ??????
  beta  = 1.0; // ??????
  
  // =================================================================
  
  /*
    Aw = (LU)w = d

    a_1 b_1  0   0   0      d_1
    c_2 a_2 b_2  0   0      d_2
     0  c_3 a_3 b_3  0   =  d_3
     0   0  c_4 a_4 b_4     d_4
     0   0   0  c_5 a_5     d_5
  */
  
  // Passos 1-3
  double a[N + 1], b[N + 1], c[N + 1], d[N + 1];
  
  double x;
  
  // Passos 4-9
  double l[N + 1], u[N + 1], z[N + 1], w[N + 2];

  // Específico
  double K[N + 2];
  
  // ============================| Começo
  time_req   = clock();
  inter_time = clock();

  // Passo 1
  h = (B - A)/(N + 1);
  
  // Passo 1.1 >> Vetor de Ks
  K[0]   = K_function(A);
  K[N+1] = K_function(B);

  for(int i = 1; i <= N; i++)
    K[i] = K_function(A + (i - 0.5) * h);

  a[1] =   (K_half(K[0], K[1]) + K_half(K[1], K[2])) / (h*h);
  b[1] = - K_half(K[1], K[2]) / (h*h);
  
  if(debug_const == 1)
    cout << "Passo 1: " << (double)(clock() - inter_time)/CLOCKS_PER_SEC << endl;
  inter_time = clock();

  // Passo 2
  for(int i = 2; i < N; i++)
  {
    x = A + (i - 0.5) * h;
    
    c[i] = - K_half(K[i-1], K[i]) / (h*h);
    a[i] =  (K_half(K[i-1], K[i]) + K_half(K[i], K[i+1])) / (h*h);
    b[i] = - K_half(K[i], K[i+1]) / (h*h);

    d[i] = q(x);
  }
  
  // Passo 3
  x = A + (N - 0.5) * h;
  
  c[N] = -  K_half(K[N-1], K[N]) / (h*h);
  a[N] =   (K_half(K[N-1], K[N]) + K_half(K[N], K[N+1])) / (h*h);

  // Fronteira com Neumann
  bool changeFronteira = true;
  if(changeFronteira)
  {
    d[1] = q(x) + (K_half(K[0], K[1]) * alpha)/(h*h);
    d[N] = q(x) + (K_half(K[N], K[N+1]) * beta) / (h*h);
  } else {
    
  }


  if(debug_const == 1)
    cout << "Passo 3: " << (double)(clock() - inter_time)/CLOCKS_PER_SEC << endl;
  inter_time = clock();

  // Passo 4
  l[1] = a[1];
  u[1] = b[1] / a[1];
  z[1] = d[1] / l[1];
  
  // Passo 5
  for(int i = 2; i < N; i++)
  {
    l[i] = a[i] - (c[i] * u[i-1]);
    u[i] = b[i] / l[i];
    z[i] = (d[i] - (c[i] * z[i-1])) / l[i];
  }
  
  // Passo 6
  l[N] = a[N] - (c[N] * u[N-1]);
  z[N] = (d[N] - (c[N] * z[N-1])) / l[N];
  
  if(debug_const == 1)
    cout << "Passo 6: " << (double)(clock() - inter_time)/CLOCKS_PER_SEC << endl;

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
  string name_file = "1D_finVolMet/1D_finVolMet_" + to_string(N) + ".txt";

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
}

// Main
int main()
{
  cout << "Versão SEM paralelização" << endl;
  cout << fixed << setprecision(12);
  // for(int i = 100; i < 1000; i+=100)
  //   FinVol(i*100);

  FinVol(25);

  return 0;
}

// g++ -o 1D_EDP-CF 1D_EDP-CF.cpp && ./1D_EDP-CF
