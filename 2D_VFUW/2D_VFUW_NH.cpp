/* ================================ | Código Exercício Chemetov-Henrique | ===================================== */
//
// Mescla a implementação do método de volumes finitos para equações elípticas com adição do método UpWind para
// resolver equações hiperbólicas de um poluente. Aqui, tenta-se a implementação de um reservatório com um esque-
// ma do tipo a quarter of the five spot.
// 
// Tenha em mente um reservatório com fluxo governado pelo sistema
// 
// ∇ . u = q em Ω    
//     p = g em ∂Ω_p 
// u . ñ = z em ∂Ω_u 
// 
// com pressão relacionada à vel. de Darcy u = (- K ∇ p) em Ω.
// 
// Então, dada a velocidade u, resolver o transporte de um 'contaminante'
// 
// ∂_t + ∇ . (uC) = 0        em Ω
// C(x,0)         = C_0(x)   em Ω
// C(x,t)         = C_D(x,t) em ∂Ω⁻ = {x ∈ ∂Ω | u.ñ_∂Ω < 0}
//
// Onde C é a concentração do poluente, C_0 a condição inicial e C_D é a concentração nas fronteiras de entrada.
// 
// 
//
/* ============================================================================================================= */

// Bibliotecas
#include <Eigen/Dense>
#include <algorithm>
#include <iostream>
#include <unistd.h>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <fftw3.h>
#include <random>
#include <string>
#include <cmath>
#include <ctime>

using namespace std;
using namespace Eigen;

const double PI  = 3.141592653589793238463;    //value of pi
const double PI2 = 9.869604401089358618834;    //value of pi^2

class GaussianField {
  private:
    int size;
    double beta;
    MatrixXd field;
    VectorXd x_grid, y_grid;

  public:
    GaussianField(int n, double b, unsigned int seed = 1) : size(n), beta(b) {
      // Initialize grids
      x_grid = VectorXd::LinSpaced(size, 0, 1);
      y_grid = VectorXd::LinSpaced(size, 0, 1);
      
      // Generate random field
      generate_field(seed);
    }

  void generate_field(unsigned int seed) {
    // Initialize random number generators
    mt19937 gen(seed);
    normal_distribution<double> dist(0.0, 1.0);

    // Allocate FFTW arrays
    fftw_complex *in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * size * size);
    fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * size * size);
    
    // Create FFTW plan
    fftw_plan plan = fftw_plan_dft_2d(size, size, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

    // Generate spectral coefficients
    for (int i = 0; i < size; ++i) {
      for (int j = 0; j < size; ++j) {
        // Compute wave numbers
        double kx = 2 * M_PI * (i <= size/2 ? i : i - size)/size;
        double ky = 2 * M_PI * (j <= size/2 ? j : j - size)/size;
        double k = sqrt(kx*kx + ky*ky);

        // Handle k=0 case
        if (i == 0 && j == 0) {
            in[i*size + j][0] = 0.0;
            in[i*size + j][1] = 0.0;
            continue;
        }

        // Power spectrum
        double power = pow(k, -(beta + 2));

        // Generate complex Gaussian noise
        double real_part = dist(gen) * sqrt(power/2);
        double imag_part = dist(gen) * sqrt(power/2);
        
        in[i*size + j][0] = real_part;
        in[i*size + j][1] = imag_part;
      }
    }

    // Perform inverse FFT
    fftw_execute(plan);

    // Store real part of result
    field.resize(size, size);
    for (int i = 0; i < size; ++i) {
      for (int j = 0; j < size; ++j) {
          field(i,j) = out[i*size + j][0];
      }
    }

    // Normalize
    double std_dev = sqrt((field.array() - field.mean()).square().sum() / (size*size));
    field /= std_dev;

    // Clean up
    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);
  }

  double interpolate(double x, double y) {
    // Find grid indices
    int i = static_cast<int>(x * (size-1));
    int j = static_cast<int>(y * (size-1));
    
    // Bilinear interpolation weights
    double dx = x * (size-1) - i;
    double dy = y * (size-1) - j;
    
    // Boundary checks
    i = min(max(i, 0), size-2);
    j = min(max(j, 0), size-2);
    
    // Bilinear interpolation
    return (1-dx)*(1-dy)*field(i,j) + 
           dx*(1-dy)*field(i+1,j) + 
           (1-dx)*dy*field(i,j+1) + 
           dx*dy*field(i+1,j+1);
  }

  const MatrixXd& get_field() const { return field; }
  const VectorXd& get_x_grid() const { return x_grid; }
  const VectorXd& get_y_grid() const { return y_grid; }
};

double K_function(GaussianField gf, double x, double y) {
  double csi = gf.interpolate(x, y);
  return exp(4.5*csi);
}

double q(double x, double y) {
  if(x == 0.02 && y == 0.02)
    return  1.0;
  else if(x == 0.98 && y == 0.98)
    return -1.0;
  else
    return  0.0;
}

// Funciona para ambos 'i' e 'j'
double K_half(double K_1, double K_2) {
  return 2.0 * (K_1 * K_2) / (K_1 + K_2);
}

// Função exata, se houver
double exact_solution(double x, double y) {
  return 0.0;
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

  // Passo 1.1 >> Geração do campo gaussiano
  GaussianField gf(M, .5);

  if(debug_const == 1)
    cout << "Passo 1: " << (double)(clock() - inter_time)/CLOCKS_PER_SEC << endl;
  inter_time = clock();

  // Passo 2 >> Montagem da Matriz
  for(int i = 1; i <= M; i++){
    for (int j = 1; j <= N; j++){
      k = i + (j - 1)*M;

      // Passo 2.1 >> Diagonais exteriores (y-direction)
      diag[0][k] = (j > 1) ? -K_half(K_function(gf,i*hx,j*hy),K_function(gf,i*hx,(j+1)*hy)) / hy2 : 0.0;
      diag[4][k] = (j < N) ? -K_half(K_function(gf,i*hx,(j-1)*hy),K_function(gf,i*hx,j*hy)) / hy2 : 0.0; 

      // Passo 2.2 >> Diagonais internas (x-direction)
      diag[1][k] = (i > 1) ? -K_half(K_function(gf,i*hx,j*hy),K_function(gf,(i+1)*hx,j*hy)) / hx2 : 0.0;
      diag[3][k] = (i < M) ? -K_half(K_function(gf,(i-1)*hx,j*hy),K_function(gf,i*hx,j*hy)) / hx2 : 0.0;

      // Passo 2.3 >> Diagonal Central
      diag[2][k] = - diag[0][k] - diag[1][k] - diag[3][k] - diag[4][k];
    }
  }

  // Passo 2.4 >> Introduzindo solução de Dirichlet local em (?,?)
  int k_ref = 1 + (1 - 1)*M;

  diag[2][k_ref] = 1.0;                  // Diagonal principal igual a 1
  diag[4][k_ref] = diag[3][k_ref] = 0.0; // Resto zerado
  diag[1][k_ref] = diag[0][k_ref] = 0.0;

  d[k_ref] = 200.0;

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

  // Passo 5.4 >> Fazer uma matriz de pressao

  double** p = new double*[M+1];
  for(int i = 0; i <= M; i++)
    p[i] = new double[N+1];

  for(int i = 1; i <= M; i++){
    for(int j = 1; j<= N; j++){
      k = i + (j - 1)*M;
      p[i][j] = w[k];
    }
  }


  if(debug_const == 1)
    cout << "Passo 5: " << (double)(clock() - inter_time)/CLOCKS_PER_SEC << endl;
  inter_time = clock();

  // Passo 6 >> Campo de Velocidades
  double*** Vel = new double**[M+1];
  for(int i = 0; i <= M; i++){
    Vel[i] = new double*[N+1];
    for(int j = 0; j <= N; j++){
      Vel[i][j] = new double[4];

      // Passo 7.2 >> Inicialização com 0.0
      for (int k = 0; k < 4; k++)
        Vel[i][j][k] = 0.0;
    }
  }

  double*** campoVel = new double**[M+1];
  for(int i = 0; i <= M; i++){
    campoVel[i] = new double*[N+1];
    for(int j = 0; j <= N; j++){
      campoVel[i][j] = new double[2];

      // Passo 7.2 >> Inicialização com 0.0
      for (int k = 0; k < 2; k++)
        campoVel[i][j][k] = 0.0;
    }
  }

  double maxVel = 0.0;
  
  for(int i = 1; i <= M; i++){
    for(int j = 1; j <= N; j++){
      k = i + (j - 1)*M;

      double Khy, Khx;
      double abs_vel;

      // Passo 6.1 >> Norte + Sul
      Khx = K_half(K_function(gf,(i-1)*hx,j*hy),K_function(gf,i*hx,j*hy));
      Vel[i][j][3] = (i > 1) ? -(Khx * (p[i][j] - p[i-1][j])) / hx : 0.0; // S

      maxVel = (fabs(Vel[i][j][3]) > maxVel) ? fabs(Vel[i][j][3]) : maxVel;

      Khx = K_half(K_function(gf,i*hx,j*hy),K_function(gf,(i+1)*hx,j*hy));
      Vel[i][j][2] = (i < M) ? -(Khx * (p[i+1][j] - p[i][j])) / hx : 0.0; // N

      maxVel = (fabs(Vel[i][j][2]) > maxVel) ? fabs(Vel[i][j][2]) : maxVel;

      // Média da velocidade em y
      campoVel[i][j][0] = (Vel[i][j][2] + Vel[i][j][3]) / 2.0;

      // Passo 6.2 >> Leste + Oeste
      Khy = K_half(K_function(gf,i*hx,j*hy),K_function(gf,i*hx,(j+1)*hy));
      Vel[i][j][0] = (j < N) ? -(Khy * (p[i][j+1] - p[i][j])) / hy : 0.0; // L

      maxVel = (fabs(Vel[i][j][0]) > maxVel) ? fabs(Vel[i][j][0]) : maxVel;

      Khy = K_half(K_function(gf,i*hx,(j-1)*hy),K_function(gf,i*hx,j*hy));
      Vel[i][j][1] = (j > 1) ? -(Khy * (p[i][j] - p[i][j-1])) / hy : 0.0; // O

      maxVel = (fabs(Vel[i][j][1]) > maxVel) ? fabs(Vel[i][j][1]) : maxVel;
      
      // Média da velocidade em x
      campoVel[i][j][1] = (Vel[i][j][0] + Vel[i][j][1]) / 2.0;

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
  double CFL = .5;
  double dt  = CFL * min(hx, hy) / maxVel;

  // Passo 7.1 >> Criar matriz de concentração Co
  int T = 20;

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

  if(debug_const == 1)
    cout << "Passo 7: " << (double)(clock() - inter_time)/CLOCKS_PER_SEC << endl;
  inter_time = clock();

  // Passo 8 >> Modelo upwind #
  // Co[n+1,i,j] = Co[n,i,j] - (Δt/Δx)(F⁺ - F⁻) - (Δt/Δy)(G⁺ - G⁻)

  for(int t = 1; t <= T; t++){
    
    // Passo 8.1 >> Pré-modelo
    double Kh;
    double Velx[M+2];
    double Vely[N+2];

    for(int i = 1; i <= M; i++){
      for(int j = 1; j <= N; j++){
        Co[t][i][j] = Co[t-1][i][j];
      }
    }

    // Eixo x
    for(int j = 1; j <= N; j++){
      // Numa linha em i
      Velx[0] = 0.0;
      Velx[M] = 0.0;

      for(int i = 1; i < M; i++){
        Kh = K_half(K_function(gf,i*hx,j*hy),K_function(gf,(i+1)*hx,j*hy));
        Velx[i] = -Kh * (p[i+1][j] - p[i][j]) / hx;
      }


      for(int i = 1; i <= M; i++){
        Co[t][i][j] -= (dt/hx)*Co[t-1][i][j]*(
          max(Velx[i],0.0) - max(Velx[i-1],0.0) +
          min(Velx[i],0.0) - min(Velx[i-1],0.0));

        Co[t][i][j] -= (dt/hx)*(max(Velx[i-1],0.0)*(Co[t-1][i][j] - Co[t-1][i-1][j])
                  + min(Velx[i],0.0)*(Co[t-1][i+1][j] - Co[t-1][i][j]));
      }
    }

    // Eixo y
    for(int i = 1; i <= M; i++){
      // Numa linha em j
      Vely[0] = 0.0;
      Vely[N] = 0.0;

      for(int j = 1; j < N; j++){
        Kh = K_half(K_function(gf,i*hx,j*hy),K_function(gf,i*hx,(j+1)*hy));
        Vely[j] = -Kh * (p[i][j+1] - p[i][j]) / hy;
      }

      for(int j = 1; j <= N; j++){
        Co[t][i][j] -= (dt/hy)*Co[t-1][i][j]*(
          max(Vely[j],0.0) - max(Vely[j-1],0.0) +
          min(Vely[j],0.0) - min(Vely[j-1],0.0));

        Co[t][i][j] -= (dt/hy)*(max(Vely[j-1],0.0)*(Co[t-1][i][j] - Co[t-1][i][j-1])
                  + min(Vely[j],0.0)*(Co[t-1][i][j+1] - Co[t-1][i][j]));
      }
    }

    cout << Velx[1] << endl;

    // Injeção
    Co[t][xinj][yinj]   = 1.0;
    Co[t][xinj-1][yinj] = 1.0;
    Co[t][xinj][yinj-1] = 1.0;

    // Co[t][xinj][ysuc]   = 1.0;
    // Co[t][xsuc][yinj]   = 1.0;

    // Extração
    Co[t][xsuc+1][ysuc] = Co[t][xsuc][ysuc];
    Co[t][xsuc][ysuc+1] = Co[t][xsuc][ysuc];
    
    // Passo 8.4 >> Entre 1 e 0, meu parceiro
    //Co[t][i][j] = fmax(0.0, fmin(1.0, Co[t][i][j]));
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
  string name_file = "2D_VFUW_db/2D_VF_" + to_string(DMR) + ".txt";

  file.open(name_file);

  file << "Tabela de valores para implementação com " << DMR << " subintervalos" << endl;
  file << fixed << setprecision(12);
  file << "Tempo de execução = " << Tempo_TOTAL << " seg.\n\n";
  file << "x;y;f(x,y);Vx;Vy" << endl;

  for (int i = 1; i <= DMR; i++)
    somaErro += abs(w[i] - valor_real[i]);

  file << fixed << setprecision(12);
  for(int i = 1; i <= M; i++){
    for(int j = 1; j <= N; j++){
      k = i + (j - 1)*M;
      file << Al + (i - 0.5) * hx  << ";";
      file << Ga + (j - 0.5) * hy  << ";";
      file << p[i][j]              << ";";
      file << campoVel[i][j][1]    << ";";
      file << campoVel[i][j][0]    << endl;
    }
  }

  file.close();

  // Impressão do campo de permeabilidades
  name_file = "2D_VFUW_db/KField_" + to_string(DMR) + ".txt";

  file.open(name_file);

  file << "Tabela de valores do campo de permeabilidade " << DMR << " subintervalos" << endl;
  file << fixed << setprecision(12);
  file << "x;y;K(x,y)" << endl;


  for (int i = 1; i <= DMR; i++)
    somaErro += abs(w[i] - valor_real[i]);

  file << fixed << setprecision(12);
  for(int i = 1; i <= M; i++){
    for(int j = 1; j <= N; j++){
      file << Al + (i - 0.5) * hx  << ";";
      file << Ga + (j - 0.5) * hy  << ";";
      file << K_function(gf, Al + (i - 0.5)*hx, Ga + (j - 0.5)*hy) << endl;
    }
  }

  file.close();

  // Binario dos timestamps
  name_file = "2D_VFUW_db/2D_UW_" + to_string(DMR) + ".txt";
  file.open(name_file);  

  file << "Tabela de Concentracao com " << DMR << " subintervalos" << endl;
  file << fixed << setprecision(12);
  file << "t;x;y;Co(x,y)" << endl;

  for(int t = 0; t <= T; t++){
    for(int i = 1; i <= M; i++){
      for(int j = 1; j <= N; j++){
        file << t*dt                 << ";";
        file << Al + (i - 0.5) * hx  << ";";
        file << Ga + (j - 0.5) * hy  << ";";
        file << Co[t][i][j]          << endl;
      }
    }
  }
  
  file.close();

  cout << "(" << DMR << ") Erro médio: " << somaErro / DMR << " em " << Tempo_TOTAL << "s" << endl;

  // ============================| Cleanup

  // Cleanup diag array
  for(int i = 0; i < 5; i++)
      delete[] diag[i];
  delete[] diag;

  // Cleanup campoVel
  for (int i = 0; i < 2; i++) {
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

// g++ -std=c++11 -O3 -I/usr/include/eigen3 -I/usr/include/ 2D_VFUW_NH.cpp -lfftw3 -lm -o 2D_VFUW_NH && ./2D_VFUW_NH