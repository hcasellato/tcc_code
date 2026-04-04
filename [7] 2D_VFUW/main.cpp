/* ===================================== | Código Final Chemetov-Henrique | ===================================== */
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
/* ============================================================================================================== */

#include <string>

const double PI  = 3.141592653589793238463; //value of pi
const double PI2 = 9.869604401089358618834; //value of pi^2

#include "permeabilidade.cpp"
#include "VFHiperbolica.cpp"
#include "VFEliptica.cpp"

int main(int argc, char const *argv[])
{
  int M = 50;
  int N = 50;
  int type_k = 0;

  ofstream file;
  file.open("2D_VFUW_db/debug_perm.txt");

  file << fixed << setprecision(12);
  file << "x;y;K(x,y)" << endl;
  for(int i = 1; i <= M; i++){
    for(int j = 1; j <= N; j++){
      file << (i) * (1.0/M)  << ";";
      file << (j) * (1.0/N)  << ";";
      file << K_function((i)*(1.0/M), (j)*(1.0/N), type_k) << endl;
    }
  }

  file.close();
  
  double** p = VFE(M, N, M, "2D_VFUW_db/debug_VFE.txt", type_k);
  
  VFH(p, M, N, "2D_VFUW_db/debug_VFH.txt", type_k);
  
  // Limpeza da memória alocada por VFE
  for(int i = 0; i <= M; i++) {
      delete[] p[i];
  }
  delete[] p;
  
  return 0;
}

// g++ main.cpp -lm -o main && ./main