#include <R.h>
#include "nrutil.h"

#define FLIP(a,b) {temp=(a);(a)=(b);(b)=temp;}

// a[1..n][1..n] is the input matrix. 
// b[1..n][1..m] is input containing the m right-hand side vectors. 
// On// output, a is replaced by its matrix inverse, and b is replaced by the corresponding set 
// of solution vectors.

void SolveHReg(double a[3][3],double b[3][1])

{
	int i, icol, irow, j, k, l, ll, p, m;
	double Amax, dtemp, PivotInverse, temp;
	int ixc[3];
	int ixr[3];
	int IP[3];
	p=3;
	m=1;
	irow = 0;
	icol = 0;
	
	for (j=0;j<p;j++) IP[j]=0;
	for (i=0;i<p;i++) {
		Amax=0.0;
		for (j=0;j<p;j++)
			if (IP[j] != 1)
				for (k=0;k<p;k++) {
					if (IP[k] == 0) {
						if (fabs(a[j][k]) >= Amax) {
							Amax=fabs(a[j][k]);
							irow=j;
							icol=k;
						}
					}
				}
		++(IP[icol]);
		if (irow != icol) {
			for (l=0;l<p;l++) FLIP(a[irow][l],a[icol][l])
			for (l=0;l<m;l++) FLIP(b[irow][l],b[icol][l])
		}
		ixr[i]=irow;
		ixc[i]=icol;
// can not occur, a[icol][icol] == 0.0, for this problem
		PivotInverse=1.0/a[icol][icol];
		a[icol][icol]=1.0;
		for (l=0;l<p;l++) a[icol][l] *= PivotInverse;
		for (l=0;l<m;l++) b[icol][l] *= PivotInverse;
		for (ll=0;ll<p;ll++)
			if (ll != icol) {
				dtemp=a[ll][icol];
				a[ll][icol]=0.0;
				for (l=0;l<p;l++) a[ll][l] -= a[icol][l]*dtemp;
				for (l=0;l<m;l++) b[ll][l] -= b[icol][l]*dtemp;
			}
	}
	for (l=p-1;l>=0;l--) {
		if (ixr[l] != ixc[l])
			for (k=0;k<p;k++)
				FLIP(a[k][ixr[l]],a[k][ixc[l]]);
	}
}
#undef FLIP
//#undef NRANSI


