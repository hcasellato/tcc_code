/* ================================ | Código Exercício Chemetov-Henrique | ===================================== */
//
// Código de Exercício, baseado nos apontamentos da Prof. Vanessa sobre o problema de PVC para ODEs, assim como
// no livro de Análise Numérica do Burden (Cap 11).
// Usa a implementação de Diferenças finitas linear (pag 769 - Burden)
//
/* ============================================================================================================= */

// Bibliotecas
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <ctime>

using namespace std;

// Function of interest:

// (0) | f(y,y',x) = y" = 4y - 4x , 0 <= x <= 1
//     | y(0) = 0
//     | y(1) = 2

// (1) | f(y,y',x) = y" = (1/2)*y - (1/2)*e^(-0.2x)

// (2) | f(y,y',x) = y" = (-2/x)y' + (2/x^2)y + sin(ln(x))/x^2, 1 <= x <= 2
//     | y(1) = 1
//     | y(2) = 2

// =================================================================
// ESCREVA AQUI AS FUNÇÕES AUXILIARES DA FORMA:
// f = y" = p(x)y' + q(x)y + r(x)

// Atualmente no problema 2
double p_function(double x)
{
  return -2.0/x;
}

double q_function(double x)
{
  return 2.0/(pow(x,2));
}

double r_function(double x)
{
  return (sin(log(x)))/(pow(x,2));
}

// Função exata, se houver
double exact_solution(double x) {
    //double e = exp(1); 
    //return pow(e, 2)*(pow(e, 2 * x) - pow(e, -2 * x)) / (pow(e, 4) - 1) + x;
    
    double c_2 = 1.0/70 * (8 - 12*sin(log(2)) - 4*cos(log(2)));
    double c_1 = 11.0/10 - c_2;
    return c_1*x + c_2/(pow(x,2)) - 3.0/10 * (sin(log(x))) - 1.0/10 * (cos(log(x)));
    
    //return 0;
    
}

// =================================================================

// Main
int main()
{
  // ================= | Variáveis!
  double A, B;         // Limites dos Intervalos e.g. x \in [A.B]
  double alpha, beta;  // Soluções de contorno
  double h;
  double valor_real;
  
  int N; // Número de Subintervalos
 
  clock_t time_req; // Para medir o tempo de execução!
  
  // =================================================================
  // COLOQUE AQUI OS VALORES DE INTERESSE
  
  N = 19;    // Número de Subintervalos
  A = 1.0;
  B = 2.0;
  
  alpha = 1.0;
  beta  = 2.0;
  
  // =================================================================
  
  /*

    b_1 c_1  0   0  
    a_2 b_2 c_2  0  
     0  a_3 b_3 c_3 
     0   0  a_4 b_4
     
  */
  
  // Passos 1-3
  double a[N + 1];
  double b[N + 1];
  double c[N + 1];
  double d[N + 1];
  
  double x, t;
  
  // Passos 4-9
  double l[N + 1];
  double u[N + 1];
  double z[N + 1];
  double w[N + 2];

  
  // ============================| Começo
  time_req = clock();
  
  // Passo 1
  h = (B - A)/(N + 1);
  x = A + h;
  
  a[1] = 2.0 + h * h * q_function(x);
  b[1] = (h / 2.0) * p_function(x) - 1;
  
  d[1] = (-1) * h * h * r_function(x) + (1 + (h/2) * p_function(x)) * alpha;
  
  // Passo 2
  for(int i = 2; i < N; i++)
  {
    x = A + i * h;
    
    a[i] = 2.0 + h * h * q_function(x);
    b[i] =        (h / 2.0) * p_function(x) - 1;
    c[i] = (-1) * (h / 2.0) * p_function(x) - 1;

    d[i] = (-1) * h * h * r_function(x);
  }
  
  // Passo 3
  x = B - h;
  
  a[N] = 2.0 + h * h * q_function(x);
  b[N] =        (h / 2.0) * p_function(x) - 1;
  c[N] = (-1) * (h / 2.0) * p_function(x) - 1;

  d[N] = (-1) * h * h * r_function(x) + (1 - (h/2.0) * p_function(x)) * beta;
  
  
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
  
  // Passo 7
  w[0]   = alpha;
  w[N+1] = beta;
  w[N]   = z[N];

  // Passo 8
  for(int i = N - 1; i >= 1; i--)
    w[i] = z[i] - (u[i]*w[i+1]);
  

// Valor encontrado
  cout << fixed << setprecision(12);
  cout << "Tempo de execução = " << (double)(clock() - time_req)/CLOCKS_PER_SEC << " seg.\n\n";
  cout << "      x_i      |       w_i      |       y_i      |      |e_a|      " << endl;
  for (int i = 0; i <= N+1; i++)
  {
    valor_real = exact_solution(A + i*h);
    cout << A + i*h << " | " << w[i] << " | " << valor_real << " | " << abs(w[i] - valor_real) << endl;
  }
  return 0;
}
// g++ -o PVC_CF-TDM PVC_CF-TDM.cpp && ./PVC_CF-TDM
