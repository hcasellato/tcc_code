#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <string>
#include <cmath>
#include <ctime>

using namespace std;


int main()
{


  // 1: .23.56.89 :O
  // 3: 12.45.78. :L
  int DMR = 9, M = 3, N = 3;

  int final1[DMR + 1], final2[DMR + 1];

  int k = 0, g = 1, i = 0;
  
  int j = 1;
  while(j <= DMR)
  {
    final1[j] = 0;
    final2[j-1] = 0;
    j += 3;
  }

  final2[DMR] = 0;

  // while(k <= DMR)
  // {
  //   final1[k]           = k;
  //   final2[g % (DMR+1)] = g;

  //   i++;
  //   g += (i+1) % 2 + 1;
  //   k += i % 2 + 1;
  // }


  while(k <= DMR)
  {
    final1[k]           = k;

    i++;
    k += i % 2 + 1;
  }

  i=0;
  while(g <= DMR)
  {
    final2[g] = g;

    i++;
    g += (i+1) % 2 + 1;
  }

  for(int p = 0; p <= DMR; p++)
    cout << p << ": " << final1[p] << " " << final2[p] << endl;

  // int DMR = 9, M = 3;

  // string final1 = "Xxxxxxxxxx";
  // string final2 = "Xxxxxxxxxx";

  // for(int i = 1; i <= DMR - M; i++) // Sul
  //   cout << i;

  // for(int i = DMR; i > DMR - M; i--) // Contorno SUL
  //   cout << "0";
  // cout << endl;
  
  // for(int i = 1; i <= M; i++) // Contorno Norte
  //   cout << "0";

  // for(int i = M + 1; i <= DMR; i++) // Norte
  //   cout << i;
  // cout << endl;

  // for(int i = 1; i <= M; i++) // Contorno
  // {
  //   final1.replace(DMR - i + 1, 1, "0");
  //   final2.replace(i, 1, "0");
  // }

  // for(int i = 1; i <= DMR - M; i++) // Diag
  // {
  //   final1.replace(i, 1, to_string(i)); // S
  //   final2.replace(M + i, 1, to_string(M + i)); // N
  // }
  // cout << final1 << endl << final2 << endl;

  int x = 1, y = 0;
  for(int r = 1; r <= DMR; r++)
  {
    cout << x << "," << y << endl;
    x = r%M + 1;
    y = (r)%3;
  }


  return 0;
}
