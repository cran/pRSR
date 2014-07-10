#include "HRegR.h"
                                
// #define pi 3.141593

// The function RegSS() computes the regression SS of the sinusoidal model
// y = mu + A*cos(2*pi*lam) + B*sin(2*pi*lam) +eps. The inputs are the frequency
// (lam), the time points of the series (t), series length (nn), the series  
// itself(yy) and a matrix X of order n by 3.

double RegSS(double lam, double *t,  long n, double *y,  MATRIX X)
{                             
	int i,j,k, m=3;  
	double a[3][3];
	double b[3][1];
	double x[3][1];
	double SSReg=0.0;
// Here's the design matrix X
	for (i=0;i<n;i++){
	 	X[i][0] = 1.0;
	 	X[i][1] = cos(TWOPI*(lam)*t[i]);
	 	X[i][2] = sin(TWOPI*(lam)*t[i]);
	}   	
// X`X	
		 for(i=0;i<m;i++){
             for(j=0;j<m;j++){
                a[i][j]=0.0;
				 for (k=0;k<n;k++){
				 a[i][j] +=(X[k][i]*X[k][j]);
				 }
              }
		}
//  X`y
		for (i=0;i<m;i++){
			b[i][0]=0.0;
				for(j=0;j<n;j++){
				b[i][0] +=(X[j][i]*y[j]);
              }
		}
		for (i=0;i<m;i++) x[i][0]=b[i][0];
       SolveHReg(a,x); 
	  //SS Regression b`X`y
		for (i=0;i<m;i++)
			SSReg +=x[i][0]*b[i][0];
		return SSReg;

}


