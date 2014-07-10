#ifndef _HREG_H_
#define _HREG_H_

#include <R.h>
#include "nrutil.h"

#define	TWOPI 6.283185307179586

void GetHReg(double *y, double *t, double *theta);

double RegSS( double, double * ,  long n, double *, MATRIX);
void SolveHReg(double [3][3], double [3][1]);


#endif /* _HREG_H_ */


