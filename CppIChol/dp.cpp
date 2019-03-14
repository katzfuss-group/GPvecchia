#define ARMA_DONT_PRINT_ERRORS

//#include <armadillo>
#include <iostream>
#include <math.h>       /* sqrt */

using namespace std;




double dot_prod(int l1, int u1, int l2, int u2, int *row_inds, double *cells){

  double result = 0.0;

  while(l1<=u1 && l2<=u2){

    if(row_inds[l1]==row_inds[l2]) {
      result += cells[l1]*cells[l2];
      l1++; l2++;
    }
    else if(l1<l2)
      l1++;
    else
      l2++;	
  }
  
  return result;
  
}




void ichol(int* ptrs, int* inds, double* vals, int N){
  
  for( int i = 0; i<N; ++i ){
    
    for( int j = ptrs[i]; j<ptrs[i+1]; ++j ){
     
      int u1 = ptrs[i];
      int u2 = ptrs[inds[j]];
      
      double dp = dot_prod( u1, ptrs[i+1]-2, u2, ptrs[inds[j] + 1]-2, inds, vals );
     
      if( inds[j] < i ){
	vals[j] = (vals[j] - dp) / vals[ ptrs[inds[j] + 1] - 1 ];
      }
      else if( inds[j]==i ){	
	vals[j] = sqrt( vals[j] - dp );
      }
      else
	cout << "ERROR" << endl;
      
    }
    
  }

}








int main(int argc, const char **argv) {

  
  //const int N = 7;
  //const int nval = 28;
  
  //int ptrs[N+1] = {0,  1,  3,  6, 10, 15, 21, 28};
  //int inds[nval] = {0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0 ,1, 2, 3, 4, 5, 6};
  //#double vals[nval] = {8.32742556,  2.07734197,  7.95711034,  0.67105164, -0.83344795, 9.858319, -0.95986073, 1.0748726, 0.3851046, 9.54222086, -0.54250083, 0.33619485, 0.21765612, -0.2376277, 10.39601622, 0.09295549, -0.11545098, -0.02322486, 0.05334551, 0.03015019, 10.33131429, -0.02661599, 0.03305713, 0.09249502, -0.01527444, -0.00863292, 0.01050412, 10.47823142};

  const int N = 3;
  const int nval = 5;

  int inds[nval] = {0, 1, 2, 1, 2};
  int ptrs[N+1] = {0, 3, 4, 5};
  double vals[nval] = {1.0000000, 0.7408182, 0.8187308, 1.0000000, 1.0000000};


  
  ichol(ptrs, inds, vals, N);
  
  for(int i=0; i<nval; ++i){
    cout << vals[i] << endl;
  }
  
  return 0;

}
