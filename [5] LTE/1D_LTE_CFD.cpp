/* ================================ | Código Exercício Chemetov-Henrique | ===================================== */
//
// Tenta a implementação de Diferenças Finitas para os problemas unidimensionais do Cap. 2
// de Equações Lineares de Transporte do livro Numerical methods for conservation laws and related Eq. do Mishra.
//
// Também com vistas ao modelo de diferenças finitas para o Equação do Calor do livro Introduction to PDE, A Com-
// putational Approach.
// Problemas da forma: 
//
// Ut + a(x,t)Ux = 0
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
const double PI  = 3.141592653589793238463;    // value of pi
const double PI2 = 9.869604401089358618834;    // value of pi^2

// ========================================================================
// PROBLEMA:
//
// U[0](x) = U(x,0) = sen(2PIx) em [0,1], com U[0] = U[N-1] e U[N] = U[1]
// 
// Sol. ex. U(x,t) = sen(2PI(x - at))
//
// ESQUEMA EXPLICITO!!!
// ========================================================================

double exact(double x, double t, double a){
	return sin(2 * PI * (x - a*t));
}

int main(){
	// Velocidade a cte
	double a = 1.0;

	// Mesh points
	int N = 100;
	int M = 10;

	// time at...
	double t = 0.0;
	double x = 0.0;

	// Omega and step
	double A = 0.0;
	double B = 1.0;

	double dx = (B - A)/(N+1);
	double dt = 1.0/M;

	// r = dt/2dx
	double r = dt/(2*dx);

	// Values Ux and Ut
	/*
		x0 --- x1 --- x2 --- ... --- xN t0
		x0 --- x1 --- x2 --- ... --- xN t1
		x0 --- x1 --- x2 --- ... --- xN t2

		...

		x0 --- x1 --- x2 --- ... --- xN tM
	*/
	double** U = new double*[M+1];
	for(int i = 0; i <= M; i++)
		U[i] = new double[N+1];

	// Passo 1 >> Contorno em U0
	// U[0](x) = U(x,0) = sen(2PIx)
	for(int j = 0; j <= N; j++)
		U[0][j] = sin(2 * PI * (A + j*dx));

	// Passo 2 >> U(x,n+1) = -ar(U(x+1,n) - U(x-1,n)) + U(x,n)
	for(int i = 1; i <= M; i++){
		for(int j = 1; j < N; j++)
			U[i][j] = -a*r*(U[i-1][j+1] - U[i-1][j-1]) + U[i-1][j];
	}

	int T = 3;

	ofstream file;
  string name_file = "LTE_data/1D_LTE_CFD_" + to_string(N) + ".txt";

  file.open(name_file);

  file << "Tabela de valores para implementação em tempo " << T << endl;
  file << fixed << setprecision(12);
  file << "x;U(x);exact_solution" << endl;

  for(int i = 1; i <= N; i++){
    file << A + i * dx               << ";";
    file << U[T][i]                  << ";";
    file << exact(A + i * dx, 3, a) << endl;
  }

  file.close();

	// Deleção
	for(int i = 0; i <= M; i++)
		delete[] U[i];
	delete[] U;

	return 0;
}

// g++ -o 1D_LTE_CFD 1D_LTE_CFD.cpp && ./1D_LTE_CFD