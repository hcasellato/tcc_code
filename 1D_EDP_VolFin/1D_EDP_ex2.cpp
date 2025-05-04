/* ================================ | Código Exercício 2 Chemetov-Henrique | ===================================== */
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
#include <vector>
#include <cmath>
#include <ctime>

#include "pbPlots.hpp"
#include "supportLib.hpp"

using namespace std;

// =================================================================
// PROBLEMA:
//
// - d/dx (K dp/dx) = 0 em          x in [0,1]
// 
// p(0) = 0 e p(1) = 1
// 
// K(x) = E                         se 0 <= x <= r
//      = 1                         se r <  x <= 1
//
// e solução exata:
// 
// p(x) = x/(r - E*r + E)           se 0 <= x <= r
// p(x) = E(x-1)/(r - E*r + E) + 1  se r <  x <= 1
// 
// com E = 1/7 e r = x_estrela = 1/3
// 
// 
// =================================================================

// Não gosto de variáveis globais, mas paciência
double epsilon   = 1.0/7;
double x_estrela = 1.0/3;

double K_function(double x) {
  if(x < x_estrela)
    return epsilon;
  else
    return 1;

  return -1; // para indicar erro de condição
             // e não falhar o retorno
}

double q(double x) {
  return 0;
}

double K_half(double K_1, double K_2) {
  // Pode ser pensada em um array de Ks, mas
  // depende da futura implementação
  
  return 2 * (K_1 * K_2) / (K_1 + K_2);
}

// Função exata, se houver
double exact_solution(double x) {
  if(x < x_estrela)
    return x / (x_estrela - epsilon*x_estrela + epsilon);
  else
    return 1+(epsilon*(x-1))/(x_estrela-epsilon*x_estrela+epsilon);

  return -1; // para indicar erro de condição
             // e não falhar o retorno
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
 
  clock_t time_req, inter_time; // Para medir o tempo de execução!
  
  // =================================================================
  // COLOQUE AQUI OS VALORES DE INTERESSE
  
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

  a[1] =  (K_half(K[0], K[1]) + K_half(K[1], K[2])) / (h*h);
  b[1] = - K_half(K[1], K[2]) / (h*h);

  // Condição de Neumann
  d[1] = q(x) + (K_half(K[0], K[1]) * A)/(h*h);
  
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

  // Condição de Neumann
  d[N] = q(x) + (K_half(K[N], K[N+1]) * B) / (h*h);
  
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
  string name_file = "1D_finVolMet_ex2/1D_finVolMet_" + to_string(N) + ".txt";

  file.open(name_file);

  file << "Tabela de valores para implementação com " << N << " subintervalos" << endl;
  file << fixed << setprecision(12);
  file << "Tempo de execução = " << Tempo_TOTAL << " seg.\n\n";
  file << "x_i;w_i;y_i;|e_a|;" << endl;

  vector<double> x_vector;
  vector<double> real_vector;
  vector<double> w_vector(w, w + (N+2));

  for (int i = 0; i <= N+1; i++)
  {
    valor_real = exact_solution(A + i*h);
    file << A + i*h << ";" << w[i] << ";" << valor_real << ";" << abs(w[i] - valor_real) << ";" << endl;
    somaErro += abs(w[i] - valor_real);
    x_vector.push_back(A + i*h); // usando o loop para assign (usado no gráfico)
    real_vector.push_back(valor_real);
  }
  file.close();

  // =================================| GRÁFICO |=================================
  bool success;
  StringReference *errorMessage = CreateStringReferenceLengthValue(0, L' ');
  RGBABitmapImageReference *imageReference = CreateRGBABitmapImageReference();

  // Cria linha de valores estimados
  ScatterPlotSeries *estimado = GetDefaultScatterPlotSeriesSettings();
  estimado->xs = &x_vector;
  estimado->ys = &w_vector;
  estimado->linearInterpolation = true;
  estimado->lineType = toVector(L"solid");
  estimado->lineThickness = 1;
  estimado->color = CreateRGBAColor(0, 0, 1, 1);

  // Cria linha de valores reais
  ScatterPlotSeries *real = GetDefaultScatterPlotSeriesSettings();
  real->xs = &x_vector;
  real->ys = &real_vector;
  real->linearInterpolation = true;
  real->lineType = toVector(L"solid");
  real->lineThickness = 1;
  real->color = CreateRGBAColor(1, 0, 0, 1);

  // Cria configurações do gráfico
  ScatterPlotSettings *settings = GetDefaultScatterPlotSettings();
  settings->width = 1200;
  settings->height = 800;
  settings->autoBoundaries = true;
  settings->autoPadding = true;

  wstring wtitle = L"1D_finVolMet_" + to_wstring(N) + L"_subintervalos";
  settings->title = toVector(wtitle.c_str());
  
  //settings->title = toVector(L"");
  settings->xLabel = toVector(L"x");
  settings->yLabel = toVector(L"Pressao");
  settings->scatterPlotSeries->push_back(real);
  settings->scatterPlotSeries->push_back(estimado);

  success = DrawScatterPlotFromSettings(imageReference, settings, errorMessage);

  if(success){
    vector<double> *pngdata = ConvertToPNG(imageReference->image);
    string name_img = "1D_finVolMet_ex2_img/1D_finVolMet_img_" + to_string(N) + ".png";
    WriteToFile(pngdata, name_img);
    DeleteImage(imageReference->image);
  }else{
    cerr << "Error: ";
    for(wchar_t c : *errorMessage->string){
      wcerr << c;
    }
    cerr << endl;
  }

  FreeAllocations();

  cout << "(" << N << ") Erro médio: " << somaErro / (N+1) << " em " << Tempo_TOTAL << "s" << endl;
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

// g++ 1D_EDP_ex2.cpp pbPlots.o supportLib.o -lm -o 1D_EDP_ex2 && ./1D_EDP_ex2