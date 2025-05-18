/* ================================ | Código Exercício Chemetov-Henrique | ===================================== */
//
// Implementação de Fatoração de Crout para problemas bidimensionais usando o método de volumes finitos.
//
// Problema:
// NABLA . (- K NABLA p) = q   em OMEGA
//                    p  = p_b em OMEGA_p
//    (- K NABLA p) . n  = u_b em OMEGA_u
// 
// Solução exata: p(x,y) = cos(2πx)cos(2πy)
// Condições de Neumann homogêneas.
/* ============================================================================================================= */

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <string>
#include <cmath>
#include <ctime>

using namespace std;
const double PI  = 3.141592653589793238463;    // Valor de π
const double PI2 = 9.869604401089358618834;    // Valor de π²

// Funções do problema
double K_function(double x) { return 1.0; }

double q(double x, double y) {
    return -8.0 * PI2 * cos(2.0 * PI * x) * cos(2.0 * PI * y);
}

double K_half(double K_1, double K_2) { return 2.0 * (K_1 * K_2) / (K_1 + K_2); }

double exact_solution(double x, double y) {
    return cos(2.0 * PI * x) * cos(2.0 * PI * y);
}

// Função principal para resolver o problema
void FinVol() {
    int M = 25, N = 25, DMR = M * N;
    double A = 0.0, B = 1.0, C = 0.0, D = 1.0;
    double hx = (B - A)/M, hy = (D - C)/N;
    double hx2 = pow(hx, 2), hy2 = pow(hy, 2);

    // Alocação dinâmica para evitar estouro de pilha
    double** diag = new double*[5];
    for (int i = 0; i < 5; i++) diag[i] = new double[DMR + 1];

    double* Kx = new double[DMR + 2]; // Permeabilidade em x
    double* Ky = new double[DMR + 2]; // Permeabilidade em y

    // Inicialização de Kx e Ky
    Kx[0] = K_function(A); Kx[DMR + 1] = K_function(B);
    Ky[0] = K_function(C); Ky[DMR + 1] = K_function(D);

    for (int i = 1; i <= M; i++) {
        for (int j = 1; j <= N; j++) {
            int k = i + (j - 1) * M;
            Kx[k] = K_function(A + (i - 0.5) * hx);
            Ky[k] = K_function(C + (j - 0.5) * hy);
        }
    }

    // Montagem da matriz pentadiagonal
    for (int k = 1; k <= DMR; k++) {
        int j = (k - 1) / M + 1;
        int i_in_row = (k - 1) % M + 1;

        // Norte e Sul
        diag[0][k] = (j > 1) ? -K_half(Ky[k], Ky[k - M]) / hy2 : 0.0;
        diag[4][k] = (j < N) ? -K_half(Ky[k], Ky[k + M]) / hy2 : 0.0;

        // Oeste e Leste
        diag[1][k] = (i_in_row > 1) ? -K_half(Kx[k], Kx[k - 1]) / hx2 : 0.0;
        diag[3][k] = (i_in_row < M) ? -K_half(Kx[k], Kx[k + 1]) / hx2 : 0.0;

        // Diagonal principal
        diag[2][k] = -(diag[0][k] + diag[1][k] + diag[3][k] + diag[4][k]);
    }

    // Vetores do sistema
    double *d = new double[DMR + 1], *z = new double[DMR + 1], *w = new double[DMR + 1];
    double *valor_real = new double[DMR + 1];

    for (int i = 1; i <= M; i++) {
        for (int j = 1; j <= N; j++) {
            int k = i + (j - 1) * M;
            d[k] = q(A + (i - 0.5)*hx, C + (j - 0.5)*hy);
            valor_real[k] = exact_solution(A + (i - 0.5)*hx, C + (j - 0.5)*hy);
        }
    }

    // Fatoração de Crout
    double** LUMatrix = new double*[5];
    for (int i = 0; i < 5; i++) LUMatrix[i] = new double[DMR + 1];

    LUMatrix[2][1] = diag[2][1];
    LUMatrix[3][1] = diag[3][1] / LUMatrix[2][1];
    z[1] = d[1] / LUMatrix[2][1];

    for (int i = 2; i <= DMR; i++) {
        double sum = 0.0;
        if (i > M) sum += LUMatrix[0][i] * LUMatrix[4][i - M];
        if (i % M != 1) sum += LUMatrix[1][i] * LUMatrix[3][i - 1];
        
        LUMatrix[2][i] = diag[2][i] - sum;
        LUMatrix[3][i] = (i % M != 0) ? diag[3][i] / LUMatrix[2][i] : 0.0;
        LUMatrix[4][i] = (i <= DMR - M) ? diag[4][i] / LUMatrix[2][i] : 0.0;

        sum = d[i];
        if (i > M) sum -= LUMatrix[0][i] * z[i - M];
        if (i % M != 1) sum -= LUMatrix[1][i] * z[i - 1];
        z[i] = sum / LUMatrix[2][i];
    }

    // Substituição regressiva
    w[DMR] = z[DMR];
    for (int i = DMR - 1; i >= 1; i--) {
        w[i] = z[i];
        if (i % M != 0) w[i] -= LUMatrix[3][i] * w[i + 1];
        if (i <= DMR - M) w[i] -= LUMatrix[4][i] * w[i + M];
    }

    // Saída dos resultados
    ofstream file("2D_finVolMet/results.txt");
    file << fixed << setprecision(12);
    for (int i = 1; i <= DMR; i++) {
        int x = (i - 1) % M + 1, y = (i - 1) / M + 1;
        file << A + (x - 0.5)*hx << " " << C + (y - 0.5)*hy << " " << w[i] << endl;
    }
    file.close();

    // Liberação de memória
    for (int i = 0; i < 5; i++) {
        delete[] diag[i];
        delete[] LUMatrix[i];
    }
    delete[] diag; delete[] LUMatrix;
    delete[] Kx; delete[] Ky; delete[] d; delete[] z; delete[] w; delete[] valor_real;
}

int main() {
    FinVol();
    return 0;
}