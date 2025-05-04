/* =========================== | Código Projeto de Pesquisa Chemetov-Henrique | ================================ */

// COCLITE[1]
// Having numerical methods for (5.1) – (5.4), we can formulate an “IMPES” method to find the 
// pressure and saturation as functions of time. As is commonly believed in the reservoir simulation
// community, it is sufficient to update the total velocity less frequently than the saturation.

/* ============================================================================================================= */
#include <iostream>
#include <cmath>
#include <cstdlib>

void simple_reservoir_sim(double* saturation_matrix,
                          int*    q,                 // The source term q accounts for injection of water and
                                                     // production of oil. Numerically, we model this term as a
                                                     // sum of “delta-functions” located at the relevant wells.
                          double  T,
                          double  time_steps         // Time steps on forLoop with k, uses N for some reason, so
                                                     // I decided that would be more clear put in separate vars.
                          int     N,                 // Rows    1 ≤ i ≤ N
                          int     M,                 // Columns 1 ≤ j ≤ M
                          double  e_cal,             // positive parameter (it's mu)
                          double  a,                 // domain limit on x
                          double  b)                 // domain limit on y

{

  double relative_permeability_w[N][M]; // matrix of water's relative permeability
  double relative_permeability_o[N][M]; // matrix of oil's relative permeability
  double total_permeability[N][M];      // sum of relative permeabilities
  
  double Lambda_right[N - 1][M]         // to solve 5.9 we need Lambda_i+1/2
  double Lambda_down[N][M - 1]          // to solve 5.9 we need Lambda_j+1/2

  double dt = T / N;

  // Making the initial guess with the LAMBDA values on the relative permeability matrix. Also giving the 
  // boundary conditions before running the actual simulation! 
  for(int row = 0; row < N; row++)
  {
    for(int column = 0; column < M; column++)
    {
      relative_permeability_w[row][column] = lambda_function(saturation_matrix[row][column], 0);
      relative_permeability_o[row][column] = lambda_function(saturation_matrix[row][column], 1);
    }
  }

  for(int k = 0; k < time_steps; k++) // iterating on time steps!!! using N for some reason unknown for me now
  {
    // ------------ solve (5.8) to get Λwij and Λoij
    update_rel_perm(relative_permeability_w, e_cal, N, M, a, b);
    update_rel_perm(relative_permeability_o, e_cal, N, M, a, b);
    
    // ------------ solve (5.9) to get vxi+1/2,j and vyi,j+1/2            !!!!!!!!!!! (for loops kinda messy)
    // first we need to find the total permeability
    for(int row = 0; row < N; row++)
    {
      for(int column = 0; column < M; column++)
      {
        total_permeability[row][column] = relative_permeability_w[row][column] + relative_permeability_o[row][column];
      }
    }
    
    // now update the permeabilities on i + 1/2 (cell 0,0 represents 1/2,0 and
    //                                           cell 1.0 represents 3/2,0. 
    for(int row = 0; row < N - 1; row++)
    {
      for(int column = 0; column < M; column++)
      {
        Lambda_right[row][column] = 2*(total_permeability[row][column] * total_permeability[row+1][column])/
                                      (total_permeability[row][column] + total_permeability[row+1][column]);
      }
    }
    
    // now update the permeabilities on j + 1/2
    for(int row = 0; row < N; row++)
    {
      for(int column = 0; column < M - 1; column++)
      {
        Lambda_down[row][column] = 2*(total_permeability[row][column] * total_permeability[row][column+1])/
                                     (total_permeability[row][column] + total_permeability[row][column+1]);
      }
    }
    
    
  }
  

}

/* ========================================== | Lambda Function | ============================================== */
// Lambda functions (water and oil) on page 16 (577) COCLITE[1]
// choose 0 for water and 1 to oil
double lambda_function(double saturation, int choice)
{
  return saturation*saturation*(1 - choice) + (1 - saturation)*(1 - saturation)*(choice);
}

/* ==================================== | Update Relative Permeability | ======================================= */
// Update relative permeability matrix
void update_rel_perm(double* relative_permeability, double e_cal, int N, int M, double a, double b)
{
  double new_value, fin_diff_x, fin_diff_y;

  double h_x = a / N;
  double h_y = b / M;

  // making the update on the values inside the border
  for(i = 1; i < N - 1; i++)
  {
    for(j = 1; j < M - 1; j++)
    {
      
      fin_diff_x = (relative_permeability[i + 1][j] - 2*relative_permeability[i][j] + relative_permeability[i - 1][j])/h_x;
      
      fin_diff_y = (relative_permeability[i][j + 1] - 2*relative_permeability[i][j] + relative_permeability[i][j - 1])/h_y;
                
      relative_permeability[i][j] = relative_permeability[i][j] - e_cal(fin_diff_x + fin_diff_y);
    }
  }
      
}

/* ========================================== | Update Velocity | ============================================== */
void update_velocity()
{








}
/* =========================================== | Main Function | =============================================== */
int main()
{
  
  // Size of the matrix N x M we are dealing with
  int N = 100;
  int M = 100;
  
  double saturation_matrix[N][M];
  int    q[N][M];

  // Making saturation_matrix and q with values = 0 for all i,j
  // i don't trust C's ability to initialize all values as 0
  for(int row = 0; row < N; row++)
  {
    for(int column = 0; column < M; column++)
    {
      saturation_matrix[row][column] = 0;
      
      q[row][column] = 0; // also updating values on q to optimize
    }
  }
  
  // Exceptions on q matrix:
  q[0][0]     =  1;
  q[N-1][M-1] = -1;
  
}
