/* ================================ | Código Exercício Chemetov-Henrique | ===================================== */
//
// Código de Exercício, baseado nos exercícios do livro de Análise Num. do Burden (pg 809)
// Implementação do Método de Diferenças Finitas usando Gauss-Seidel para equações elípticas (de Poisson) 
/*

    u_xx(x,y) + u_yy(x,y) = f(x,y)   a <= x <= b ; c <= y <= d

    u(x,y) = g(x,y)                  x = a ou x = b    e    c <= y <= d
    u(x,y) = g(x,y)                  y = c ou y = d    e    a <= x <= b

*/
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

// ===================================================================
// Funçoes Aux

double f(double x, double y)
{
  return x * exp(y);
}

// Condições de contorno
double g(double x, double y)
{
  // Precisa mesmo desses ifs esquisitos?
  // tá muito feio pra ser correto
  if(x == 0)
    return 0;
  if(x == 2)
    return 2 * exp(y);

  if(y == 0)
    return x;
  if(y == 1)
    return x * exp(1);


  return -1000; // Para indicar erro
}

double u(double x, double y)
{
  return x * exp(y);
}


// ===================================================================

// Main
int main()
{
  // =================================================================
  // COLOQUE AQUI OS VALORES DE ENTRADA
  
  double a, b, c, d;   // Limites dos Intervalos

  // Intervalo x
  a = 0;
  b = 2;
  
  // Intervalo y
  c = 0;
  d = 1;

  // Tamanho da matriz
  int m, n;
  m = 5;               // >= 3
  n = 6;               // >= 3
  
  // Auxiliares
  double TOL = .0000000001; // Tolerância
  int    N   = 100;       // Nro Max de iterações

  clock_t time_req; // Para medir o tempo de execução!

  // =================================================================
  // Variáveis extras internas
  
  double h, k; // Discretização dos intervalos
  
  // Malha de pontos
  double x[n];
  double y[m];
  
  // Malha de resultado
  double w[n][m];
  
  // ???
  double lambda, mu, l;
  
  // Importante nas iterações de Gauss-Seidel
  double z, NORMA;
  
  
  // ============================| Começo
  
  // Passo 1
  h = (b - a)/n;
  k = (d - c)/m;
  
  // Passo 2 e 3
  for(int i = 1; i < n; i++)
    x[i] = a + i*h;
  
  for(int j = 1; j < m; j++)
    y[j] = a + j*k;
  
  // Passo 4
  for(int i = 1; i < n; i++)
  {
    for(int j = 1; j < m; j++)
      w[i][j] = 0;  
  }
  
  // Passo 5
  lambda = pow(h,2) / pow(k,2);
  mu     = 2 * (1 + lambda);
  l      = 1;
  
  // Passo 6
  while(l <= N) // Passos 7 a 20 implementam iterações de Gauss-Seidel
  {
    // Passo 7
    z = (-h*h*f(x[1],y[m-1]) + g(a,y[m-1]) + lambda*g(x[1],d) + lambda*w[1][m-2] + w[2][m-1])/mu;
    
    NORMA = abs(z - w[1][m-1]);
    
    w[1][m-1] = z;
    
    // Passo 8
    for(int i = 2; i <= n-2; i++)
    {
      z = (-h*h*f(x[i],y[m-1]) + lambda*g(x[i],d) + w[i-1][m-1] + w[i+1][m-1] + lambda*w[i][m-2])/mu;
      
      if(abs(z - w[i][m-1]) > NORMA)
        NORMA = abs(z - w[i][m-1]);
        
      w[i][m-1] = z;
    }
    
    // Passo 9
    z = (-h*h*f(x[n-1],y[m-1]) + g(b,y[m-1]) + lambda*g(x[n-1],d) + lambda*w[n-1][m-2] + w[n-2][m-1])/mu;
    
    if(abs(z - w[n-1][m-1]) > NORMA)
      NORMA = abs(z - w[n-1][m-1]);
    
    w[n-1][m-1] = z;
    
    // Passo 10
    for(int j = m - 2; j >= 2; j--)
    {
      // Passo 11
      z = (-h*h*f(x[1],y[j]) + g(a,y[j]) + lambda*w[1][j-1] + lambda*w[1][j+1] + w[2][j])/mu;
      
      if(abs(z - w[1][j]) > NORMA)
        NORMA = abs(z - w[1][j]);
      
      w[1][j] = z;

      // Passo 12
      for(int i = 2; i <= n-2; i++)
      {
        z = (-h*h*f(x[i],y[j]) + w[i-1][j] + lambda*w[i][j+1] + w[i+1][j] + lambda*w[i][j-1])/mu;
        
        if(abs(z - w[i][j]) > NORMA)
          NORMA = abs(z - w[i][j]);
          
        w[i][j] = z;
      }
      
      // Passo 13
      z = (-h*h*f(x[n-1],y[j]) + g(b,y[j]) + w[n-2][j] + lambda*w[n-1][j+1]+lambda*w[n-1][j-1])/mu;
      
      if(abs(z - w[n-1][j]) > NORMA)
        NORMA = abs(z - w[n-1][j]);
      
      w[n-1][j] = z;
      
    }
    
    // Passo 14
    z = (-h*h*f(x[1],y[1]) + g(a,y[1]) + lambda*g(x[1],c) + lambda*w[1][2] + w[2][1])/mu;
    
    if(abs(z - w[1][1]) > NORMA)
      NORMA = abs(z - w[1][1]);
    
    w[1][1] = z;
    
    // Passo 15
    for(int i = 2; i <= n-2; i++)
    {
      z = (-h*h*f(x[i],y[1]) + lambda*g(x[i],c) + w[i-1][1] + w[i+1][1] + lambda*w[i][2])/mu;
      
      if(abs(z - w[i][1]) > NORMA)
        NORMA = abs(z - w[i][1]);
        
      w[i][1] = z;
    }

    // Passo 16
    z = (-h*h*f(x[n-1],y[1]) + g(b,y[1]) + lambda*g(x[n-1],c) + lambda*w[n-1][2] + w[n-2][1])/mu;
    
    if(abs(z - w[n-1][1]) > NORMA)
      NORMA = abs(z - w[n-1][1]);
    
    w[n-1][1] = z;
    
    // Passo 17
    if(NORMA <= TOL)
    {
      double Tempo_TOTAL = (double)(clock() - time_req)/CLOCKS_PER_SEC;
      // Impressão dos Valores
      // Por favor, crie a pasta antes!
      ofstream file;
      string name_file = "2D_Burden/2D_EDP_CN_" + to_string((int)log10(1/TOL)) + ".txt";

      file.open(name_file);

      file << "Tabela de valores de w para implementação com tolerância de " << TOL << endl;
      file << fixed << setprecision(12);
      file << "Tempo de execução = " << Tempo_TOTAL << " seg.\n\n";

      // Passo 18
      for(int i = 1; i < n; i++)
      {
        for(int j = 1; j < m; j++)
        {
          file << w[i][j] << ";";
        }
        file << endl;
      }
      file.close();

      cout << l << " iterações com tolerância de " << TOL << " em " << Tempo_TOTAL << "s" << endl;
    
      // Passo 19
      return 0;
    }
    // Passo 20
    l++;

  }
  
  // Passo 21
  cout << "Processo mal sucedido!!!" << endl;
  cout << "(número máximo de iterações excedido)" << endl;

  return 0;
}
// g++ -o 2D_EDP_DFGS 2D_EDP_DFGS.cpp && ./2D_EDP_DFGS
