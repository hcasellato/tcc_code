/* ================================ | Código Exercício Chemetov-Henrique | ===================================== */
//
// Tenta resolver o exercício 1 dos projetos computacionais 6.8 do livro de Métodos de Volumes Finitos (Prof. Dr.
// Fabricio, pg. 91). 
//
// - ∇ . u = q      em Ω
//       p = p_b    em ∂Ω_p
//   u . n = u_b    em ∂Ω_u com u = (- K ∇ p) em Ω.
//
// ∂_t + ∇ . (uC) = 0        em Ω
// C(x,0)         = C_0(x)   em Ω
// C(x,t)         = C_D(x,t) em ∂Ω⁻ = {x ∈ ∂Ω | u.ñ_∂Ω < 0}
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
// - d/dx (K dp/dx) = -25cos(25x) em  Ω = [0,1]
//                p = x           em ∂Ω
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
  return 2 * (K_1 * K_2) / (K_1 + K_2);
}

// Função exata, se houver
double exact_solution(double x) {
  return x;
}

double co(double x) {
	return exp(-20.0*(x-2.0)*(x-2.0)) + exp(-(x-5.0)*(x-5.0));
}

// =================================================================
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
  
  alpha = 0.0;
  beta  = 1.0;
  
  // =================================================================
  
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
  
  x = A + (0.5 * h);
	d[1] = q(x) + (K_half(K[0], K[1]) * alpha)/(h*h);
  
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

	d[N] = q(x) + (K_half(K[N], K[N+1]) * beta) / (h*h);

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

  if(debug_const == 1)
    cout << "Passo 8: " << (double)(clock() - inter_time)/CLOCKS_PER_SEC << endl;

  // Passo 9 >> Velocidades
  // u = k/μ . (p[i+1] - p[i])/h
  double Vel[N + 2];
  double VelFac[N + 2];

  double face_sul = -K_half(K[0], K[1]) * (w[1] - w[0]) / h;
  double face_nor;

  double maxVel = 0.0;

  VelFac[0] = face_sul;
  Vel[0]    = face_sul; // Contorno sul

  for(int i = 1; i <= N; i++){
  	face_nor  = -K_half(K[i], K[i+1]) * (w[i+1] - w[i])/h;
  	Vel[i]    = (face_nor + face_sul) * .5;

  	VelFac[i] = face_nor;

  	face_sul  = face_nor;

    maxVel    = (abs(VelFac[i]) > maxVel) ? abs(VelFac[i]) : maxVel;
  }

  VelFac[N+1] = face_nor;
  Vel[N+1]    = face_nor; // Contorno norte

  maxVel    = (abs(VelFac[N+1]) > maxVel) ? abs(VelFac[N+1]) : maxVel;

	// Passo 10 >> UpWind
	double CFL = .9;
	double dt  = CFL * h / maxVel;

	int T = 100;

	double** Co = new double*[T+1];
	for(int i = 0; i <= T; i++)
		Co[i] = new double[N+2];

	for(int i = 0; i <= N; i++){
		x = A + (i + .5) * h;
		//Co[0][i] = 0.0;

    Co[0][i] = co(3*x+1.5);
	}

	for(int i = 1; i <= T; i++){

    //Co[i][0] = 1.0 - (double)i / T;
    Co[i][N+1]   = abs(.1 * sin(10.0*PI*i*dt) + .2/(2.0*i));
    Co[i][0] = .0; 
		for(int j = 1; j <= N; j++){

			Co[i][j]  = Co[i-1][j];
      
      Co[i][j] -= (dt/h)*Co[i-1][j]*(
        max(VelFac[j],0.0) - max(VelFac[j-1],0.0) +
        min(VelFac[j],0.0) - min(VelFac[j-1],0.0));

      Co[i][j] -= (dt/h)*(max(VelFac[j-1],0.0)*(Co[i-1][j] - Co[i-1][j-1])
        + min(VelFac[j],0.0)*(Co[i-1][j+1]-Co[i-1][j]));
		}
	}

  double Tempo_TOTAL = (double)(clock() - time_req)/CLOCKS_PER_SEC;
  double somaErro = 0;
  
  // Impressão dos Valores
  // Por favor, crie a pasta antes!
  ofstream file;
  string name_file = "1D_VFUW_dta/1D_VF_" + to_string(N) + ".txt";

  file.open(name_file);

  file << "Tabela de valores para implementação com " << N << " subintervalos" << endl;
  file << fixed << setprecision(12);
  file << "Tempo de execução = " << Tempo_TOTAL << " seg.\n\n";
  file << "x;w;y;erro;Vel" << endl;
  for (int i = 0; i <= N+1; i++)
  {
    valor_real = exact_solution(A + i*h);
    file
    << A + i*h                << ";"
    << w[i]                   << ";"
    << valor_real             << ";"
    << abs(w[i] - valor_real) << ";"
    << Vel[i]
    << endl;
    somaErro += abs(w[i] - valor_real);
  }
  file.close();

  name_file = "1D_VFUW_dta/1D_UW_" + to_string(N) + ".txt";

  file.open(name_file);

  file << "Tabela de valores para implementação com " << N << " subintervalos" << endl;
  file << fixed << setprecision(12);
  file << "Tempo de execução = " << Tempo_TOTAL << " seg.\n\n";
  file << "t;x;Co" << endl;

  for(int i = 0; i <= T; i++){
  	for(int j = 0; j <= N+1; j++){
	    file
	    << i*dt     << ";"
	    << A + j*h  << ";"
	    << Co[i][j]
	    << endl;
  	}
  }
  file.close();

  for (int i = 0; i <= T; i++)
    delete[] Co[i];
  delete[] Co;


  cout << "(" << N << ") Erro médio: " << somaErro / (N+1) << " em " << Tempo_TOTAL << "s" << endl;
}

// Main
int main()
{
  cout << "Versão SEM paralelização" << endl;
  cout << fixed << setprecision(12);

  FinVol(100);

  return 0;
}

// g++ -o 1D_VFUW 1D_VFUW.cpp && ./1D_VFUW