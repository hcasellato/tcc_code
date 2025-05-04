#include <iostream>
#include <cmath>
#include <iomanip>
#include <ctime>

using namespace std;

double K_function(double x) {
    return 2.0 + sin(25.0 * x);
}

double q(double x) {
    return -25.0 * cos(25.0 * x);
}

double K_half(double K_1, double K_2) {
    return 2.0 * (K_1 * K_2) / (K_1 + K_2);
}

double exact_solution(double x) {
    return x;
}

int main() {
    const int N = 19;
    const double A = 0.0, B = 1.0;
    const double alpha = 0.0, beta = 1.0;
    const double h = (B - A) / (N + 1);
    
    double K[N + 2];
    K[0] = K_function(A); // Left boundary
    K[N + 1] = K_function(B); // Right boundary
    
    // Corrected: Cell centers at (i - 0.5)*h
    for (int i = 1; i <= N; i++) {
        K[i] = K_function(A + (i - 0.5) * h);
    }
    
    double a[N + 1], b[N + 1], c[N + 1], d[N + 1];
    
    // Matrix Assembly
    for (int i = 1; i <= N; i++) {
        const double x = A + (i - 0.5) * h;
        
        // Interface conductivities
        const double K_left = (i == 1) ? K_half(K[0], K[1]) : K_half(K[i - 1], K[i]);
        const double K_right = (i == N) ? K_half(K[N], K[N + 1]) : K_half(K[i], K[i + 1]);
        
        a[i] = (K_left + K_right) / (h * h);
        d[i] = q(x);
        
        if (i < N) b[i] = -K_right / (h * h);
        if (i > 1) c[i] = -K_left / (h * h);
        
        // Boundary terms
        if (i == 1) d[i] += (K_left * alpha) / (h * h);
        if (i == N) d[i] += (K_right * beta) / (h * h);
    }
    
    // Crout Factorization
    double l[N + 1], u[N + 1], z[N + 1], w[N + 2];
    
    l[1] = a[1];
    u[1] = b[1] / l[1];
    z[1] = d[1] / l[1];
    
    for (int i = 2; i <= N - 1; i++) {
        l[i] = a[i] - c[i] * u[i - 1];
        u[i] = b[i] / l[i];
        z[i] = (d[i] - c[i] * z[i - 1]) / l[i];
    }
    
    l[N] = a[N] - c[N] * u[N - 1];
    z[N] = (d[N] - c[N] * z[N - 1]) / l[N];
    
    // Back substitution
    w[0] = alpha;
    w[N + 1] = beta;
    w[N] = z[N];
    
    for (int i = N - 1; i >= 1; i--) {
        w[i] = z[i] - u[i] * w[i + 1];
    }
    
    // Output
    double somaErro = 0;
    cout << fixed << setprecision(12);
    cout << "      x_i      |       w_i      |       y_i      |      |e_a|      " << endl;
    for (int i = 0; i <= N + 1; i++) {
        const double x = A + i * h;
        const double exato = exact_solution(x);
        const double erro = abs(w[i] - exato);
        cout << x << " | " << w[i] << " | " << exato << " | " << erro << endl;
        somaErro += erro;
    }
    cout << "Erro médio: " << somaErro / (N + 2) << endl;
    
    return 0;
}