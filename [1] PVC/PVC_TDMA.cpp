/* ================================ | Código Exercício Chemetov-Henrique | ===================================== */
//
// Código de Exercício, baseado nos apontamentos da Prof. Vanessa sobre o problema de PVC para ODEs.
// Usa a implementação do algoritmo de Thomas, para matrizes tridiagonais!
//
/* ============================================================================================================= */

// Bibliotecas
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>

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
  // Matriz da forma Aw = v
  
  // ================= | Variáveis!
  double a, b;         // Limites dos Intervalos e.g. x \in [a.b]
  double sol_a, sol_b; // Soluções de contorno
  double h;
  double valor_real;
  
  int N;
  
  // =================================================================
  // COLOQUE AQUI OS VALORES DE INTERESSE
  
  N = 10; // Número de Subintervalos
  a = 1;
  b = 2;
  
  sol_a = 1;
  sol_b = 2;
  
  // =================================================================
  
  double malha_x[N+1]; // Malha de x_i
  double malha_w[N+1]; // Malha de soluções

  // Para o algo. de Thomas, são necessários 3 vetores, a, b, c
  // representando as diagonais abaixo, do meio e acima tal como:
  /*
  
      b_1 c_1  0   0  
      a_2 b_2 c_2  0  
       0  a_3 b_3 c_3 
       0   0  a_4 b_4
       
  */
  // com a_1 = c_n = 0.
  
  double diag_a[N-1], diag_b[N-1], diag_c[N-1];
  double malha_v[N+1];
  
  double diag_c_new[N-1], malha_v_new[N+1]; // valores novos para Alg. de Thomas
  
  // ================= | Passo 1:
  h = (b - a)/N;
  
  malha_x[0] = a;
  malha_x[N] = b;
  
  for(int i = 1; i < N; i++)
    malha_x[i] = a + i*h;
  
  // ================= | Passo 2? e 3:
  malha_w[0] = sol_a;
  malha_w[N] = sol_b;
  
  // ================= | Passo 4:
  
  // Conseguir as diagonais a, b, c
  
  for(int i = 0; i < N-1; i++)
  {
    diag_a[i] = (-1) * (h / 2.0) * p_function(malha_x[i + 1]) - 1;

    diag_b[i] = 2.0 + h * h * q_function(malha_x[i + 1]);

    diag_c[i] = (h / 2.0) * p_function(malha_x[i + 1]) - 1;
  }

  // valores para adequação ao método (início de a e final de c)
  diag_a[0] = 0;
  diag_c[N-2] = 0;

  for(int i = 1; i < N; i++)
    malha_v[i] = -1.0*pow(h,2)*r_function(malha_x[i]);

  malha_v[1]   += (1 + (h/2.0)*p_function(malha_x[0]))*sol_a;
  malha_v[N-1] += (1 - (h/2.0)*p_function(malha_x[N]))*sol_b;
  
  // Resolver pelo Algoritmo de Thomas (também suporta paralelização)
  // Matriz da forma Aw = v
                                           
  diag_c_new[0] = diag_c[0] / diag_b[0];
  for(int i = 1; i < N-1; i++)
    diag_c_new[i] = diag_c[i] / (diag_b[i] - diag_c_new[i-1] * diag_a[i]);
  
  // e.g. para N = 4:
  // a =     a_0 a_1 a_2
  // v = v_0 v_1 v_2 v_3 v_4
  
  malha_v_new[1] = malha_v[1] / diag_b[0];
  for(int i = 2; i < N; i++)
    malha_v_new[i] = (malha_v[i] - malha_v_new[i-1] * diag_a[i-1]) / (diag_b[i-1] - diag_c_new[i-2] * diag_a[i-1]);

  // Solução:
  malha_w[N-1] = malha_v_new[N-1];
  for(int i = N-2; i > 0; i--)
    malha_w[i] = malha_v_new[i] - diag_c_new[i-1] * malha_w[i+1];
  
  for(int i = 0; i < N+1; i++)
    cout << i << " " << malha_w[i] << endl;

  // Impressão para solução
  cout << fixed << setprecision(8);
  cout << "    x_i    |    w_i     |    y_i     |    |e_a|    " << endl;
  for (int i = 0; i <= N; i++)
  {
    valor_real = exact_solution(malha_x[i]);
    cout << malha_x[i] << " | " << malha_w[i] << " | " << valor_real << " | " << abs(malha_w[i] - valor_real) << endl;
  }
  
  return 0;
}
// g++ -o PVC_TDMA PVC_TDMA.cpp && ./PVC_TDMA
