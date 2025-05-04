/* ================================ | Código Exercício Chemetov-Henrique | ===================================== */
//
// Tenta a implementação de Fatoração de Crout de forma paralela para os problemas unidimensionais do livro de
// Volumes Finitos.
// 
//
// Problemas da forma: 
//
// - d/dx (K dp/dx) = q
//
/* ============================================================================================================= */

// Bibliotecas
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <memory>
#include <string>
#include <thread>
#include <vector>
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

// Função de cálculo para K_half paralelizavel
// provavelmente terei de fazer uma matriz cabd ?
void paralell_KHalf(double** cabd, double* K, int inicio, int fim, double A, double h)
{
  double x;
  for(int i = inicio; i <= fim; i++)
  {
    x = A + (i - 0.5) * h;
    
    cabd[0][i] = - K_half(K[i-1], K[i]) / (h*h);
    cabd[1][i] =  (K_half(K[i-1], K[i]) + K_half(K[i], K[i+1])) / (h*h);
    cabd[2][i] = - K_half(K[i], K[i+1]) / (h*h);

    cabd[3][i] = q(x);
  }
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
  
  // Exclusivo para multiThreading
  int num_threads = 4;
  thread threads[num_threads];

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
  double** cabd = new double*[4];
  for (int i = 0; i < 4; i++)
    cabd[i] = new double[N+1];
  
  double x;
  
  // Passos 4-9
  double l[N + 1], u[N + 1], z[N + 1], w[N + 2];

  // Específico
  double* K = new double[N + 2];
  
  // ============================| Começo
  //cout << N << endl;
  time_req = clock();
  inter_time = clock();

  // Passo 1
  h = (B - A)/(N + 1);
  
  // Passo 1.1 >> Vetor de Ks
  K[0]   = K_function(A);
  K[N+1] = K_function(B);

  if(debug_const == 1)
    cout << "Passo 1: " << (double)(clock() - inter_time)/CLOCKS_PER_SEC << endl;
  inter_time = clock();

  for(int i = 1; i <= N; i++)
    K[i] = K_function(A + (i - 0.5) * h);

  cabd[1][1] =   (K_half(K[0], K[1]) + K_half(K[1], K[2])) / (h*h);
  cabd[2][1] = - K_half(K[1], K[2]) / (h*h);

  cabd[3][1] = q(x) + (K_half(K[0], K[1]) * alpha)/(h*h);
  
  // Passo 2
  
  threads[0] = thread(paralell_KHalf, cabd, K, 2, 1*(int)(N/num_threads), A, h);
  for(int i = 1; i < num_threads - 1; i++)
    threads[i] = thread(paralell_KHalf, cabd, K, i*(int)(N/num_threads) + 1, (i+1)*(int)(N/num_threads), A, h);
  threads[num_threads - 1] = thread(paralell_KHalf, cabd, K, (num_threads - 1)*(int)(N/num_threads) + 1, N - 1, A, h);

  for (auto& th : threads) th.join();
  

  // Passo 3
  x = A + (N - 0.5) * h;
  
  cabd[0][N] = -  K_half(K[N-1], K[N]) / (h*h);
  cabd[1][N] =   (K_half(K[N-1], K[N]) + K_half(K[N], K[N+1])) / (h*h);

  cabd[3][N] = q(x) + (K_half(K[N], K[N+1]) * beta) / (h*h);
  
  if(debug_const == 1)
    cout << "Passo 3: " << (double)(clock() - inter_time)/CLOCKS_PER_SEC << endl;
  inter_time = clock();
  // Numa análise preliminar, os passos 1-3 incluso podem ser
  // paralelizados.

  // Passo 4
  l[1] = cabd[1][1];
  u[1] = cabd[2][1] / cabd[1][1];
  z[1] = cabd[3][1] / l[1];
  
  // Passo 5
  for(int i = 2; i < N; i++)
  {
    l[i] = cabd[1][i] - (cabd[0][i] * u[i-1]);
    u[i] = cabd[2][i] / l[i];
    z[i] = (cabd[3][i] - (cabd[0][i] * z[i-1])) / l[i];
  }
  // Passo 6
  l[N] = cabd[1][N] - (cabd[0][N] * u[N-1]);
  z[N] = (cabd[3][N] - (cabd[0][N] * z[N-1])) / l[N];
  
  if(debug_const == 1)
    cout << "Passo 6: " << (double)(clock() - inter_time)/CLOCKS_PER_SEC << endl;
  inter_time = clock();

  // Passo 7
  w[0]   = alpha;
  w[N+1] = beta;
  w[N]   = z[N];

  // Passo 8
  for(int i = N - 1; i >= 1; i--)
    w[i] = z[i] - (u[i]*w[i+1]);

  double Tempo_TOTAL = (double)(clock() - time_req)/CLOCKS_PER_SEC;

  //cout << "Passo 8: " << (double)(clock() - time_req)/CLOCKS_PER_SEC << endl;
  double somaErro = 0;
  
  // Impressão dos Valores
  // Por favor, crie a pasta antes!
  ofstream file;
  string name_file = "1D_finVolMet_Pa/1D_finVolMet_Pa" + to_string(N) + ".txt";

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

  for (int i = 0; i < 4; i++)
    delete[] cabd[i];
  delete[] cabd;
  delete[] K;
}

// Main
int main()
{
  cout << "Versão com 8 threads" << endl;
  cout << fixed << setprecision(12);
  for(int i = 100; i < 1000; i+=100)
    FinVol(i*100);

  return 0;
}

// g++ -std=c++11 -pthread -o 1D_EDP-CF_Pa 1D_EDP-CF_Pa.cpp && ./1D_EDP-CF_Pa
