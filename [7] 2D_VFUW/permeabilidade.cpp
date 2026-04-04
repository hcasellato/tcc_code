/* ================================ | Código Exercício Chemetov-Henrique | ===================================== */
//
// Cria um campo de permeabilidade para ser usado nos outros códigos.
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
/* ============================================================================================================= */

// Bibliotecas
#include <algorithm>
#include <iostream>
#include <unistd.h>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <fftw3.h>
#include <random>
#include <vector>
#include <string>
#include <cmath>
#include <ctime>

using namespace std;

double K_function(double x, double y, int type) {
	double d;

	switch(type){
		case 0: // Campo homogêneo
			return 1.0;
			break;

		case 1: // Campo nH -> faixas verticais
			if(x >= .0 && x <= 1.0){
				if       (y > .2 && y <= .4){
					return 0.3;
				} else if(y > .4 && y <= .6){
					return 0.7;
				} else if(y > .6 && y <= .8){
					return 0.3;
				}
			}
			return 0.1;
			break;

		case 2: // Campo nH -> diagonal
			return exp(-fabs(x-y)-0.3);
			break;

		case 3:
			d = (x-.5)*(x-.5) + (y-.5)*(y-.5);
			return max(exp(.8*sqrt(d)+0.1)-1.1,0.1);
			break;

		default:
			return 0.0;
			break;
	}
	return 1.0;
}

// Funciona para ambos 'i' e 'j'
double K_half(double K_1, double K_2) {
  return 2 * (K_1 * K_2) / (K_1 + K_2);
}

// g++ -std=c++11 permeabilidade.cpp -o permeabilidade && ./permeabilidade.cpp