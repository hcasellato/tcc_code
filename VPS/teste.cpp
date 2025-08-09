/* ================================ | Código Exercício Chemetov-Henrique | ===================================== */
//
// Tenta mesclar a implementação de BLU Fact. para matrizes tridiagonais de blocos de problemas de volumes fini-
// tos com adição de um poluente governado pelo método upwind. Usa os códigos '2D_EDP_BLUv4' e '1D_TLE_UW'.
// 
// Tenha em mente um reservatório com fluxo governado pelo sistema
// 
// ∇ . u = q em Ω                ou = f [SOUSA, Fabricio S.]
//     p = g em ∂Ω_p
// u . ñ = z em ∂Ω_u             aqui ter atenção ao vetor ñ em "
// 
// com pressão relacionada à vel. de Darcy u = (- K ∇ p) em Ω.
// 
// Então, para u, resolver o transporte de um 'contaminante'
// 
// ∂_t + ∇ . (uC) = 0        em Ω
// C(x,0)         = C_0(x)   em Ω
// C(x,t)         = C_D(x,t) em ∂Ω⁻ = {x ∈ ∂Ω | u.ñ_∂Ω < 0}
//
// Onde C é a concentração do poluente, C_0 a condição inicial e C_D é a concentração nas fronteiras de entrada.
// 
// Nesse caso,
// 
// ∇ . u =  -8 PI^2 cos(2PI x)cos(2PI y) em Ω
//     p =  u                            em ∂Ω_p
// u . ñ =  0                            em ∂Ω_u
// 
// com permeabilidade absoluta K = 1, e, com porosidade uniforme φ = 1,
// 
// ∂_t + ∇ . (uC) = 0 em Ω
// C(x,0)         = 0 em Ω
// C(x,t)         = 1 em ∂Ω⁻ = {x ∈ ∂Ω | u.ñ_∂Ω < 0}
// 
// Que os deuses tenham misericórida do que estou prestes a digitar.
//
/* ============================================================================================================= */

// Bibliotecas
#include <algorithm>
#include <iostream>
#include <unistd.h>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <random>
#include <string>
#include <cmath>
#include <ctime>

using namespace std;
const double PI  = 3.141592653589793238463;    //value of pi
const double PI2 = 9.869604401089358618834;    //value of pi^2

double K_function(double x, double y) {
  return 1.0;
}

int test_x = 1;
int test_y = 1;

double q(double x, double y) {
  return (4.0 * test_x + 4.0 * test_y) * PI2 * cos(2.0 * PI * x * test_x) * cos(2.0 * PI * y * test_y);
}

// Funciona para ambos 'i' e 'j'
double K_half(double K_1, double K_2) {
  return 2 * (K_1 * K_2) / (K_1 + K_2);
}

// Função exata, se houver
double exact_solution(double x, double y) {
  return cos(2.0 * PI * x * test_x) * cos(2.0 * PI * y * test_y);
}

void Copy_A_to_Atil(double **A, double **Atil, int block_index){
  int TB    = 25;
  int start = (block_index - 1) * TB; // DMR - TB = 600

  Atil[1][1] = A[2][start + 1]; // 601
  Atil[1][2] = A[3][start + 1];

  for(int i = 2; i < TB; i++){
    Atil[i][i-1] = A[1][start + i];
    Atil[i][i]   = A[2][start + i];
    Atil[i][i+1] = A[3][start + i];
  }

  Atil[TB][TB-1] = A[1][start + TB]; // 625
  Atil[TB][TB]   = A[2][start + TB];

}

// Resolve sistemas lineares por decomposição LU
// Dx = f => L(Ux) = f => Lz = f => Ux = z
/*
  INPUT:
    Matriz  D[0..TB]
    Vetor   F[0..TB]
    Vetor   X[0..TB]

  O código abaixo é a fatoração LU do livro de Análi-
  se Numérica do Burden, p. 450, algoritmo 6.4.
*/
void LU_solve(double **D, double* X, double* F){
  int TB = 25;
  double sum;

  double* z = new double[TB+1];

  // Passo 0 >> Fazer matrizes L e U
  double** L = new double*[TB + 1];
  double** U = new double*[TB + 1];
  for(int i = 0; i <= TB; i++){
    L[i] = new double[TB+1];
    U[i] = new double[TB+1];
    for(int j = 1; j <= TB; j++){
      L[i][j] = 0.0;
      U[i][j] = 0.0;
    }
    U[i][i] = 1.0; // diagonal = 1
  }

  // Passo 1
  L[1][1] = D[1][1];

  for (int j = 1; j <= TB; j++) {
    // Calculate column j of L
      for (int i = j; i <= TB; i++) {
        sum = 0.0;
        for (int k = 1; k < j; k++) {
          sum += L[i][k] * U[k][j];
        }
      L[i][j] = D[i][j] - sum;
    }
    if (L[j][j] == 0.0) {
      throw runtime_error("Error: Singular matrix encountered in LU_solve. Division by zero.");
    }
    
    // Calculate row j of U
    for (int i = j + 1; i <= TB; i++) {
      sum = 0.0;
      for (int k = 1; k < j; k++) {
        sum += L[j][k] * U[k][i];
      }
      U[j][i] = (D[j][i] - sum) / L[j][j];
    }
  }

  // Passo 7.1 >> Solucionar Lz = F
  z[1] = F[1] / L[1][1];

  // Passo 7.2
  for(int i = 2; i <= TB; i++){
    sum = 0.0;
    for(int j = 1; j <= i - 1; j++)
      sum += L[i][j] * z[j];
    z[i] = (F[i] - sum)/L[i][i];
  }

  // Passo 9.1 >> Substituição regressiva
  X[TB] = z[TB];

  // Passo 9.2
  for(int i = TB - 1; i >= 1; i--){
    sum = 0.0;
    for(int j = i+1; j <= TB; j++)
      sum += U[i][j] * X[j];
    X[i] = (z[i] - sum);
  }

  // Cleanup memory
  delete[] z;
  for (int i = 0; i <= TB; i++) {
    delete[] L[i];
    delete[] U[i];
  }
  delete[] L;
  delete[] U;
}

void FinVol(int debug_const)
{
  // ================= | Variáveis!
  int M, N; // Dimensões M,N da matriz e
  int DMR;  // Dimensão da Matriz Resultante (DMR)
  int TB;   // Tamanho do bloco

  double Al, Be, Ga, De;   // (x,y) \in [Al.Be] X [Ga,De]
  double hx, hy, hx2, hy2; // Pulo entre x_i e x_{i+1}
 
  clock_t time_req, inter_time; // Para medir o tempo de execução!
  
  // =================================================================
  // COLOQUE AQUI OS VALORES DE INTERESSE
  
  Al = Ga = 0.0;
  Be = De = 1.0;
  
  M   = 25;
  N   = 25;
  TB  = 25;

  DMR = M * N;

  // ============================| DIAGONAIS

  double** diag = new double*[5];
  for(int i = 0; i < 5; i++)
    diag[i] = new double[DMR+1];

  int k; // aritmética k = i + (j - 1)M
  double x, y;
  double d[DMR+1];
  double valor_real[DMR+1];   // Vetor de valores reais da solução
  
  // RESPOSTA!
  double w[DMR + 1];
  
  // ============================| Começo 
  time_req = inter_time = clock();

  // Passo 1
  hx = (Be - Al)/M;
  hy = (De - Ga)/N;

  hx2 = hx * hx;
  hy2 = hy * hy;

  for(int i = 1; i <= M; i++){
    for(int j = 1; j <= N; j++){
      k = i + (j - 1)*M;

      // Passo 1.1 >> Vetor 'd' e das soluções exatas 
      d[k]          = q(Al + (i - 0.5) * hx, Ga + (j - 0.5) * hy);
      valor_real[k] = exact_solution(Al + (i - 0.5) * hx, Ga + (j - 0.5) * hy);
    }
  }

  if(debug_const == 1)
    cout << "Passo 1: " << (double)(clock() - inter_time)/CLOCKS_PER_SEC << endl;
  inter_time = clock();

  // Passo 2 >> Montagem da Matriz
  for(int i = 1; i <= M; i++){
    for (int j = 1; j <= N; j++){
      k = i + (j - 1)*M;

      // Passo 2.1 >> Diagonais exteriores (y-direction)
      diag[0][k] = (j > 1) ? -K_half(K_function(i*hx,j*hy),K_function(i*hx,(j+1)*hy)) / hy2 : 0.0;
      diag[4][k] = (j < N) ? -K_half(K_function(i*hx,(j-1)*hy),K_function(i*hx,j*hy)) / hy2 : 0.0; 

      // Passo 2.2 >> Diagonais internas (x-direction)
      diag[1][k] = (i > 1) ? -K_half(K_function(i*hx,j*hy),K_function((i+1)*hx,j*hy)) / hx2 : 0.0;
      diag[3][k] = (i < M) ? -K_half(K_function((i-1)*hx,j*hy),K_function(i*hx,j*hy)) / hx2 : 0.0;

      // Passo 2.3 >> Diagonal Central
      diag[2][k] = - diag[0][k] - diag[1][k] - diag[3][k] - diag[4][k];
    }
  }

  // Passo 2.4 >> Introduzindo solução de Dirichlet local em (1,1)

  int k_ref = 1;

  diag[2][k_ref] = 1.0;                  // Diagonal principal igual a 1
  diag[4][k_ref] = diag[3][k_ref] = 0.0; // Resto zerado
  
  d[k_ref] = valor_real[k_ref];

  if(debug_const == 1)
    cout << "Passo 2: " << (double)(clock() - inter_time)/CLOCKS_PER_SEC << endl;
  inter_time = clock();

  // ============================| Block LU Decomposition

  // Passo 3 >> Fazendo as matrizes
  double*** Atil = new double**[TB+1];
  double*** G    = new double**[TB+1];
  double*** B    = new double**[TB+1];

  double**  C    = new double*[TB+1];
  double**  Z    = new double*[TB+1];

  double**  b_d  = new double*[TB+1]; // bloco de d
  double**  b_w  = new double*[TB+1]; // bloco de w

  double*   b_i  = new double[TB+1]; // vetor intermediario

  // Alocando número de blocos
  for(int i = 0; i <= TB; i++){
    Atil[i] = new double*[TB+1];
    G[i]    = new double*[TB+1];
    B[i]    = new double*[TB+1];

    // Alocando linhas
    C[i]    = new double[TB+1];
    Z[i]    = new double[TB+1];
    
    b_d[i]  = new double[TB+1];
    b_w[i]  = new double[TB+1];

    // Alocando linhas
    for(int j = 0; j <= TB; j++){
      Atil[i][j] = new double[TB+1];
      G[i][j]    = new double[TB+1];
      B[i][j]    = new double[TB+1];
      
      // Inicializando colunas
      C[i][j]    = 0.0;
      Z[i][j]    = 0.0;
      
      b_d[i][j]  = 0.0;
      b_w[i][j]  = 0.0;
      
      // Inicializando colunas
      for(int k = 0; k <= TB; k++){
        Atil[i][j][k] = 0.0; 
        G[i][j][k]    = 0.0; 
        B[i][j][k]    = 0.0; 
      }
    }
  }

  // Passo 3.1.1 >> Ã1 = A1
  Copy_A_to_Atil(diag, Atil[1], 1);

  // Passo 3.1.2 >> Ã1 * G1 = C1 => G1
  for(int i = 1; i <= TB; i++){
    // Pré: fazer C1 (vou ter q fazer isso sempre)
    C[i][i] = diag[4][i];

    LU_solve(Atil[1], G[1][i], C[i]);
  }

  // Passo 3.2.1 >> Ãi = Ai - Bi * G[i-1], i = 2, 3, ..., n-1
  for(int b = 2; b < TB; b++){
    int bli = (b-1)*TB;

    // Primeiro copio os valores de diag para Atil no bloco
    Copy_A_to_Atil(diag, Atil[b], b);

    for(int i = 1; i <= TB; i++){
      for(int j = 1; j <= TB; j++){
        // Depois subtraio
        Atil[b][i][j] -= diag[0][bli+i] * G[b-1][i][j];
      }
    }

    // Passo 3.3 >> Ãi * Gi = Ci => Gi,    i = 2, 3, ..., n-1
    for(int i = 1; i <= TB; i++){
      // Pré: fazer Ci
      C[i][i] = diag[4][bli+i];

      LU_solve(Atil[b], G[b][i], C[i]);
    }
  }

  // Passo 3.2.2 >> Ãn = An - Bn * G[n-1], i = n
  int bli = (TB-1)*TB;

  Copy_A_to_Atil(diag, Atil[TB], TB);

  for(int i = 1; i <= TB; i++){
    for(int j = 1; j <= TB; j++){
      Atil[TB][i][j] -= diag[0][bli+i] * G[TB-1][i][j];
    }
  }

  if(debug_const == 1)
    cout << "Passo 3: " << (double)(clock() - inter_time)/CLOCKS_PER_SEC << endl;
  inter_time = clock();

  // Pré: fazer bloco de d
  for(int b = 1; b <= TB; b++){
    int bli = (b-1)*TB;
    for(int j = 1; j <= TB; j++)
      b_d[b][j] = d[bli + j];
  }

  // Passo 4.1 >> Ã1 * z1 = d1 => z1
  LU_solve(Atil[1], Z[1], b_d[1]);

  // Passo 4.2 >> Ãi * zi = di - Bi * z[i-1] => zi, i = 2, 3, ..., n
  for(int b = 2; b <= TB; b++){

    int bli = (b-1)*TB;
    // Fazendo di - Bi * z[i-1]
    for(int j = 1; j <= TB; j++)
      b_i[j] = b_d[b][j] - diag[0][bli+j] * Z[b-1][j];

    // Resolvendo Ãi * zi = bi = di - Bi * z[i-1]
    LU_solve(Atil[b], Z[b], b_i);
  }

  if(debug_const == 1)
    cout << "Passo 4: " << (double)(clock() - inter_time)/CLOCKS_PER_SEC << endl;
  inter_time = clock();

  // Passo 5.1 >> xn = zn
  for(int j = 1; j <= TB; j++)
    b_w[TB][j] = Z[TB][j];

  // Passo 5.2 >> xi = zi - Gi * x[i+1],     i = n-1, n-2, ..., 1
  for(int b = TB-1; b >= 1; b--){
    for(int i = 1; i <= TB; i++){
      // Calcular Gi * x[i+1]
      double sum = 0.0;
      for(int j = 1; j <= TB; j++)
        sum += G[b][i][j] * b_w[b+1][j];

      // xi = zi - Gi * x[i+1]
      b_w[b][i] = Z[b][i] - sum;
    }
  }

  // Passo 5.3 >> Transformar block_w to w:
  for(int i = 1; i <= TB; i++){
    int bli = (i-1)*TB;

    for(int j = 1; j <= TB; j++)
      w[bli+j] = b_w[i][j];
  }

  // Passo 6 >> Campo de Velocidades
  double** campoVel = new double*[4];
  for(int i = 0; i < 4; i++)
    campoVel[i] = new double[DMR+1];

  double maxVel = 0.0;
  
  for(int i = 1; i <= M; i++){
    for(int j = 1; j <= N; j++){
      k = i + (j - 1)*M;

      double Khy, Khx;
      double abs_vel;

      // Inicializando com 0
      for(int ii = 0; ii < 4; ii++)
      	campoVel[ii][k] = 0.0;

      // Passo 6.1 >> Leste + Oeste
      Khx = K_half(K_function(i*hx,j*hy),K_function((i+1)*hx,j*hy));

      abs_vel = fabs(Khx * (w[k] - w[k-1])) / hx;
      if (abs_vel > maxVel) maxVel = abs_vel;

      campoVel[0][k] += (i > 1) ? -(Khx * (w[k] - w[k-1])) / (hx * 2) : 0.0;
      campoVel[2][k] += (i > 1) ? -(Khx * (valor_real[k] - valor_real[k-1])) / (hx * 2) : 0.0;

      Khx = K_half(K_function((i-1)*hx,j*hy),K_function(i*hx,j*hy));

      abs_vel = fabs(Khx * (w[k+1] - w[k])) / hx;
      if (abs_vel > maxVel) maxVel = abs_vel;

      campoVel[0][k] += (i < M) ? -(Khx * (w[k+1] - w[k])) / (hx * 2) : 0.0;
      campoVel[2][k] += (i < M) ? -(Khx * (valor_real[k+1] - valor_real[k])) / (hx * 2) : 0.0;

      // Passo 6.2 >> Norte + Sul

      Khy = K_half(K_function(i*hx,j*hy),K_function(i*hx,(j+1)*hy));

      abs_vel = fabs(Khy * (w[k] - w[k-M])) / hy;
      if (abs_vel > maxVel) maxVel = abs_vel;

      campoVel[1][k] += (j > 1) ? -(Khy * (w[k] - w[k-M])) / (hy * 2) : 0.0;
      campoVel[3][k] += (j > 1) ? -(Khy * (valor_real[k] - valor_real[k-M])) / (hy * 2) : 0.0;

			Khy = K_half(K_function(i*hx,(j-1)*hy),K_function(i*hx,j*hy));

      abs_vel = fabs(Khy * (w[k+M] - w[k])) / hy;
      if (abs_vel > maxVel) maxVel = abs_vel;

      campoVel[1][k] += (j < N) ? -(Khy * (w[k+M] - w[k])) / (hy * 2) : 0.0;
      campoVel[3][k] += (j < N) ? -(Khy * (valor_real[k+M] - valor_real[k])) / (hy * 2) : 0.0;
  	}
  }

	// Passo 6.3 >> Calcular velocidades máximas nas faces
  if (maxVel == 0.0) {
    maxVel = 1.0;
  }

  if(debug_const == 1)
    cout << "Passo 6: " << (double)(clock() - inter_time)/CLOCKS_PER_SEC << endl;
  inter_time = clock();

  // Passo 7 >> Pré-modelo
  double dt = min(hx, hy) / maxVel;

  // Passo 7.1 >> Criar matriz de concentração Co
  int T = 50;

  double*** Co = new double**[T];
	for(int i = 0; i < T; i++){
		Co[i] = new double*[M+1];
  	for(int j = 0; j <= M; j++){
  		Co[i][j] = new double[N+1];

  		// Passo 7.2 >> Inicialização com 0.4
      for (int k = 0; k <= N; k++)
        Co[i][j][k] = 0.2;
  	}
  }

  // Passo 7.3 >> Criação de fluxo de injeção
  int xinj = 13;
  int yinj = 13;

  int xsuc = 7;
  int ysuc = 13;

  double Q_inj = 0.5;
  double vazao = Q_inj / (hx * hy);

  if(debug_const == 1)
    cout << "Passo 7: " << (double)(clock() - inter_time)/CLOCKS_PER_SEC << endl;
  inter_time = clock();

  // Passo 8 >> Modelo upwind #
  // Co[n+1,i,j] = Co[n,i,j] - (Δt/Δx)(F⁺ - F⁻) - (Δt/Δy)(G⁺ - G⁻)

  for(int t = 1; t < T; t++){
  	for(int i = 1; i <= M; i++){
  		for(int j = 1; j <= N; j++){
  			k = i + (j - 1)*M;

  			// Passo 8.1 >> Pré-modelo
				double KhxP, KhxM, KhyP, KhyM;
        double salto_plus, salto_min;

        double bordaCima, bordaBaixo;

				Co[t][i][j] = Co[t-1][i][j];

				// Passo 8.2 >> Fluxo no eixo x
				KhxP = K_half(K_function(i*hx,j*hy),K_function((i+1)*hx,j*hy));
				KhxM = K_half(K_function((i-1)*hx,j*hy),K_function(i*hx,j*hy));

				bordaBaixo = (i > 1) ? max(-(KhxP*(w[k]-w[k-1]))/hx,0.0) : 0.0;
				bordaCima  = (i < M) ? max(-(KhxM*(w[k+1]-w[k]))/hx,0.0) : 0.0;
				
				Co[t][i][j] -= Co[t-1][i][j]*(dt/hx)*(bordaBaixo - bordaCima);

				bordaBaixo = (i > 1) ? min(-(KhxP*(w[k]-w[k-1]))/hx,0.0) : 0.0;
				bordaCima  = (i < M) ? min(-(KhxM*(w[k+1]-w[k]))/hx,0.0) : 0.0;
				
				Co[t][i][j] -= Co[t-1][i][j]*(dt/hx)*(bordaBaixo - bordaCima);

				salto_plus = (i < M) ? Co[t-1][i+1][j] - Co[t-1][i][j] : 0.0;
				salto_min  = (i > 1) ? Co[t-1][i][j] - Co[t-1][i-1][j] : 0.0;

				Co[t][i][j] -= (dt/hx)*(max(-(KhxM*(w[k+1]-w[k]))/hx,0.0)*salto_min 
					+ min(-(KhxP*(w[k]-w[k-1]))/hx,0.0)*salto_plus);

				// Passo 8.2 >> Fluxo no eixo y
				KhyP = K_half(K_function(i*hx,j*hy),K_function(i*hx,(j+1)*hy));
				KhyM = K_half(K_function(i*hx,(j-1)*hy),K_function(i*hx,j*hy));

				bordaBaixo = (j > 1) ? max(-(KhyP*(w[k]-w[k-M]))/hy,0.0) : 0.0;
				bordaCima  = (j < N) ? max(-(KhyM*(w[k+M]-w[k]))/hy,0.0) : 0.0;
				
				Co[t][i][j] -= Co[t-1][i][j]*(dt/hy)*(bordaBaixo - bordaCima);

				bordaBaixo = (j > 1) ? min(-(KhyP*(w[k]-w[k-M]))/hy,0.0) : 0.0;
				bordaCima  = (j < N) ? min(-(KhyM*(w[k+M]-w[k]))/hy,0.0) : 0.0;
				
				Co[t][i][j] -= Co[t-1][i][j]*(dt/hy)*(bordaBaixo - bordaCima);

				salto_plus = (j < N) ? Co[t-1][i][j+1] - Co[t-1][i][j] : 0.0;
				salto_min  = (j > 1) ? Co[t-1][i][j] - Co[t-1][i][j-1] : 0.0;

				Co[t][i][j] -= (dt/hy)*(max(-(KhyM*(w[k+M]-w[k]))/hx,0.0)*salto_min 
					+ min(-(KhyP*(w[k]-w[k-M]))/hx,0.0)*salto_plus);

				// Passo 8.3 >> Injeção de ? contaminantes
				// int cont_times = 3;
				// for(int c = 0; c < cont_times; c++){
				// 	xinj = rand() % 25 + 1;
				// 	yinj = rand() % 25 + 1;

				// 	if(i == xinj && j == yinj){
				// 		Co[t][i][j] += dt*vazao;
				// 	}
				// }

        if(i == xinj && j == yinj){
          Co[t][i][j] += dt*vazao;
        }

        // Extração
        if(i == xsuc && j == ysuc){
          Co[t][i][j] -= dt*vazao*0.7;
        }

				// Passo 8.4 >> Entre 1 e 0, meu parceiro
				Co[t][i][j] = fmax(0.0, fmin(1.0, Co[t][i][j]));

  		}
  	}
  }

  if(debug_const == 1)
    cout << "Passo 8: " << (double)(clock() - inter_time)/CLOCKS_PER_SEC << endl;
  inter_time = clock();
  
  // ============================| Print
  double Tempo_TOTAL = (double)(clock() - time_req)/CLOCKS_PER_SEC;
  double somaErro = 0.0;
  
  // Impressão dos Valores
  // Por favor, crie a pasta antes!
  ofstream file;
  string name_file = "VPS_db/VPS_FVM_" + to_string(DMR) + ".txt";

  file.open(name_file);

  file << "Tabela de valores para implementação com " << DMR << " subintervalos" << endl;
  file << fixed << setprecision(12);
  file << "Tempo de execução = " << Tempo_TOTAL << " seg.\n\n";
  file << "x;y;f(x,y);exact_solution;difference;Vx;Vy;VxE;VyE" << endl;

  for (int i = 1; i <= DMR; i++)
    somaErro += abs(w[i] - valor_real[i]);

  file << fixed << setprecision(12);
  for(int i = 1; i <= M; i++){
    for(int j = 1; j <= N; j++){
      k = i + (j - 1)*M;
      file << Al + (i - 0.5) * hx  << ";";
      file << Ga + (j - 0.5) * hy  << ";";
      file << w[k]                 << ";";
      file << valor_real[k]        << ";";
      file << valor_real[k] - w[k] << ";";
      file << campoVel[0][k]       << ";";
      file << campoVel[1][k]       << ";";
      file << campoVel[2][k]       << ";";
      file << campoVel[3][k]       << endl;
    }
  }

  file.close();

  // Binario dos timestamps
  name_file = "VPS_db/VPS_FVM_" + to_string(DMR) + ".bin";
  ofstream fileBin(name_file, ios::binary);
  
  if (!fileBin)
    cerr << "Erro ao criar arquivo!" << endl;

  int dims[3] = {T, M, N};
  fileBin.write(reinterpret_cast<char*>(dims), sizeof(dims));
  for (int t = 0; t < T; t++) {
    for (int i = 1; i <= M; i++) {
      fileBin.write(reinterpret_cast<char*>(&Co[t][i][1]), N * sizeof(double));
    }
  }
  
  fileBin.close();

  cout << "(" << DMR << ") Erro médio: " << somaErro / DMR << " em " << Tempo_TOTAL << "s" << endl;

  // ============================| Cleanup

  // Cleanup diag array
  for(int i = 0; i < 5; i++)
      delete[] diag[i];
  delete[] diag;

  // Cleanup campoVel
  for (int i = 0; i < 4; i++) {
      delete[] campoVel[i];
  }
  delete[] campoVel;

  // Cleanup Co
  for (int t = 0; t < T; t++) {
      for (int i = 0; i <= M; i++) {
          delete[] Co[t][i];
      }
      delete[] Co[t];
  }
  delete[] Co;

  // Cleanup Atil, G, and B arrays
  for(int i = 0; i <= TB; i++) {

    // Cleanup Atil[i], G[i], B[i]
    if (Atil[i] != nullptr) {
      for(int j = 0; j <= TB; j++)
        delete[] Atil[i][j];
      delete[] Atil[i];
    }
    if (G[i] != nullptr) {
      for(int j = 0; j <= TB; j++)
        delete[] G[i][j];
      delete[] G[i];
    }
    if (B[i] != nullptr) {
      for(int j = 0; j <= TB; j++)
        delete[] B[i][j];
      delete[] B[i];
    }

    // Cleanup C[i], Z[i], b_d[i], b_w[i]
    delete[] C[i];
    delete[] Z[i];
    delete[] b_d[i];
    delete[] b_w[i];
  }

  // Delete the top-level arrays
  delete[] Atil;
  delete[] G;
  delete[] B;
  delete[] C;
  delete[] Z;
  delete[] b_d;
  delete[] b_w;
  delete[] b_i;

  if(debug_const == 1)
    cout << "Cleanup complete" << endl;

}

int main()
{
  cout << "Versão 4 da Fatoração LU em Bloco" << endl;
  cout << fixed << setprecision(12);
  FinVol(1);

  return 0;
}

// g++ teste.cpp -lm -o teste && ./teste 