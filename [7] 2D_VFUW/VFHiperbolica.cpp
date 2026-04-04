/* ================================ | Código Exercício Chemetov-Henrique | ===================================== */
//
// Implementação de Block LU factorization para matrizes em bloco tridiagonais para os problemas bidimensionais
// do livro de Volumes Finitos (Fabricio). Consultar também ANALYSIS OF NUMERICAL METHODS, Isaacson e Keller.
// 
// Essa versão é a refatoração de 2D_EDP_BLU_v4, tentando otimizar espaço e performance, além de fazer o código
// mais legível.
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
#include <unistd.h>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <string>
#include <cmath>
#include <ctime>

using namespace std;

void VFH(double** p, int M, int N, string name_file, int type_k){
  double Al, Be, Ga, De;
  double hx  = 1.0/M;
  double hy  = 1.0/N;
  int k;

  Al = Ga = 0.0;
  Be = De = 1.0;

  double vel;
  double maxVel = 0.0;
  
  for(int i = 1; i <= M; i++){
    for(int j = 1; j <= N; j++){
      k = i + (j - 1)*M;

      double Khy, Khx;

      // Passo 6.1 >> Norte + Sul
      Khx = K_half(K_function((i-1)*hx,j*hy,type_k),K_function(i*hx,j*hy,type_k));
      vel = (i > 1) ? -(Khx * (p[i][j] - p[i-1][j])) / hx : 0.0; // S

      maxVel = (fabs(vel) > maxVel) ? fabs(vel) : maxVel;

      Khx = K_half(K_function(i*hx,j*hy,type_k),K_function((i+1)*hx,j*hy,type_k));
      vel = (i < M) ? -(Khx * (p[i+1][j] - p[i][j])) / hx : 0.0; // N

      maxVel = (fabs(vel) > maxVel) ? fabs(vel) : maxVel;

      // Passo 6.2 >> Leste + Oeste
      Khy = K_half(K_function(i*hx,j*hy,type_k),K_function(i*hx,(j+1)*hy,type_k));
      vel = (j < N) ? -(Khy * (p[i][j+1] - p[i][j])) / hy : 0.0; // L

      maxVel = (fabs(vel) > maxVel) ? fabs(vel) : maxVel;

      Khy = K_half(K_function(i*hx,(j-1)*hy,type_k),K_function(i*hx,j*hy,type_k));
      vel = (j > 1) ? -(Khy * (p[i][j] - p[i][j-1])) / hy : 0.0; // O

      maxVel = (fabs(vel) > maxVel) ? fabs(vel) : maxVel;
    }
  }

  // Passo 6.3 >> Calcular velocidades máximas nas faces
  if (maxVel == 0.0) {
    maxVel = 1.0;
  }

  // Passo 7 >> Pré-modelo
  double CFL = .9;
  double dt  = CFL * min(hx, hy) / maxVel;

  // Passo 7.1 >> Criar matriz de concentração Co
  int T = 100;
  int Ttemp = 20;
  int Tcounter = 1;

  double*** Co = new double**[T+1];
	for(int i = 0; i <= T; i++){
		Co[i] = new double*[M+2];
  	for(int j = 0; j <= M+1; j++){
  		Co[i][j] = new double[N+2];

  		// Passo 7.2 >> Inicialização com 0.0
      for (int k = 0; k <= N+1; k++)
        Co[i][j][k] = 0.0;
  	}
  }

  double*** CoTemp = new double**[2];
  for(int i = 0; i < 2; i++){
    CoTemp[i] = new double*[M+2];
    for(int j = 0; j <= M+1; j++){
      CoTemp[i][j] = new double[N+2];

      // Passo 7.2 >> Inicialização com 0.0
      for (int k = 0; k <= N+1; k++)
        CoTemp[i][j][k] = 0.0;
    }
  }

  // Passo 7.3 >> Criação de fluxo de injeção
  int xinj = 1;
  int yinj = 1;

  int xsuc = M;
  int ysuc = N;

  // Injeção
  Co[0][xinj][yinj]   = 1.0;
  Co[0][xinj-1][yinj] = 1.0;
  Co[0][xinj][yinj-1] = 1.0;

  // Extração
  Co[0][xsuc+1][ysuc] = Co[0][xsuc][ysuc];
  Co[0][xsuc][ysuc+1] = Co[0][xsuc][ysuc];

  // Passo 8 >> Modelo upwind
  // Co[n+1,i,j] = Co[n,i,j] - (Δt/Δx)(F⁺ - F⁻) - (Δt/Δy)(G⁺ - G⁻)

  for(int i = 1; i <= M; i++){
    for(int j = 1; j <= N; j++){
      CoTemp[0][i][j] = Co[0][i][j];
    }
  }

  for(int t = 1; t <= T*Ttemp + 1; t++){
    // Passo 8.1 >> Pré-modelo
    double Kh;
    double Velx[M+2];
    double Vely[N+2];

    for(int i = 1; i <= M; i++){
      for(int j = 1; j <= N; j++){
        CoTemp[1][i][j] = CoTemp[0][i][j];
      }
    }

    if(t % Ttemp == 0){
      for(int i = 1; i <= M; i++){
        for(int j = 1; j <= N; j++){
          Co[Tcounter][i][j] = CoTemp[0][i][j];
        }
      }
      Tcounter++;
    }

    // Eixo x
    for(int j = 1; j <= N; j++){
      // Numa linha em i
      Velx[0] = 0.0;
      Velx[M] = 0.0;

      for(int i = 1; i < M; i++){
        Kh = K_half(K_function(i*hx,j*hy,type_k),K_function((i+1)*hx,j*hy,type_k));
        Velx[i] = -Kh * (p[i+1][j] - p[i][j]) / hx;
      }

      for(int i = 1; i <= M; i++){
        CoTemp[1][i][j] -= (dt/hx)*CoTemp[0][i][j]*(
          max(Velx[i],0.0) - max(Velx[i-1],0.0) +
          min(Velx[i],0.0) - min(Velx[i-1],0.0));

        CoTemp[1][i][j] -= (dt/hx)*(max(Velx[i-1],0.0)*(CoTemp[0][i][j] - CoTemp[0][i-1][j])
                  + min(Velx[i],0.0)*(CoTemp[0][i+1][j] - CoTemp[0][i][j]));
      }
    }

    // Eixo y
    for(int i = 1; i <= M; i++){
      // Numa linha em j
      Vely[0] = 0.0;
      Vely[N] = 0.0;

      for(int j = 1; j < N; j++){
        Kh = K_half(K_function(i*hx,j*hy,type_k),K_function(i*hx,(j+1)*hy,type_k));
        Vely[j] = -Kh * (p[i][j+1] - p[i][j]) / hy;
      }

      for(int j = 1; j <= N; j++){
        CoTemp[1][i][j] -= (dt/hy)*CoTemp[0][i][j]*(
          max(Vely[j],0.0) - max(Vely[j-1],0.0) +
          min(Vely[j],0.0) - min(Vely[j-1],0.0));

        CoTemp[1][i][j] -= (dt/hy)*(max(Vely[j-1],0.0)*(CoTemp[0][i][j] - CoTemp[0][i][j-1])
                  + min(Vely[j],0.0)*(CoTemp[0][i][j+1] - CoTemp[0][i][j]));
      }
    }

    // Injeção
    if(t < 80){      
      CoTemp[1][xinj][yinj]   = 1.0;
      CoTemp[1][xinj-1][yinj] = 1.0;
      CoTemp[1][xinj][yinj-1] = 1.0;
    }

    // Extração
    CoTemp[1][xsuc][ysuc]   = 0.0;
    CoTemp[1][xsuc+1][ysuc] = CoTemp[1][xsuc][ysuc];
    CoTemp[1][xsuc][ysuc+1] = CoTemp[1][xsuc][ysuc];
    
		// Passo 8.4
    for(int i = 1; i <= M; i++){
      for(int j = 1; j <= N; j++){
        CoTemp[0][i][j] = CoTemp[1][i][j];
      }
    }
  }

  cout << "VFH Completa!" << endl;

  ofstream file;
  file.open(name_file);  

  file << fixed << setprecision(12);
  file << "t;x;y;Co(x,y)" << endl;

  for(int t = 0; t <= T; t++){
    for(int i = 1; i <= M; i++){
      for(int j = 1; j <= N; j++){
        file << (t*Ttemp)*dt         << ";";
        file << Al + (i - 0.5) * hx  << ";";
        file << Ga + (j - 0.5) * hy  << ";";
        file << Co[t][j][i]          << endl;
      }
    }
  }

  // Cleanup Co
  for (int t = 0; t < T; t++) {
    for (int i = 0; i <= M; i++) {
      delete[] Co[t][i];
    }
    delete[] Co[t];
  }
  delete[] Co;
}