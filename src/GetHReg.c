// Fits 4 parameter harmonic regression using OLS
// y = mu + A cos(2 pi f t) + B sin(2 pi f t) + e
// Input Variables:
//  *y is pointer to a vector containing the time series (missing values excluded)
//  *t is pointer to a vector containing the time points (in full data case: 1, ..., n)
//  *theta, theta[0] = n, number of data points
//  *theta, theta[1] = nF, number of frequencies enumerated
//
// Output Variables:
//   *theta[0] contains -2LLR (the likelihood ratio statistic)
//   *theta[1] optimal frequency
//
#include "HRegR.h"
void GetHReg(double *y, double *t, double *theta)
{                             
	int i,j, n, nF; 
	VECTOR lam;
	MATRIX X;
	VECTOR Reg;
	double ymean, sumY, SSTot, SSE, MaxRegSS, lamOpt, aF, an;
	SSTot = 0.0;
	sumY = 0.0;
	n= (int) theta[0];
	nF= (int) theta[1];
	an= theta[0];
	X=Matrix(n,3);
        Reg = Vector(nF);
	lam = Vector(nF);
	aF = 2*nF+1;
	for (i=0; i<nF; i++)
		lam[i]=(i+1.0)/aF;
// regression sum of squares for each frequency  
	for(i=0;i<nF;i++)
		Reg[i]=RegSS(lam[i],t,n,y,X);  
// maximum regression sum of squares	
	MaxRegSS = Reg[0];
	lamOpt = lam[0];
	for (j=1;j<nF;j++)
		if(Reg[j]>MaxRegSS) {
			MaxRegSS = Reg[j];
			lamOpt = lam[j];
			}
// total sum of squares
	for (i = 0;i<n;i++)
		sumY +=y[i];
	ymean=sumY/n;
	for(i = 0;i<n;i++)
		SSTot +=((y[i]-ymean)*(y[i]-ymean));
// error sum of squares, remove mean	
	SSE = SSTot - (MaxRegSS - an*ymean*ymean);	
// likelihood-ratio test statistic
    theta[0] = -an*log(SSE/SSTot);
    theta[1] = lamOpt;
// free pointers
	free_matrix(X);
	free_vector(lam);
	free_vector(Reg);
	return;
}

