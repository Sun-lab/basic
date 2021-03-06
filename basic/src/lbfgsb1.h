/*
 *  lbfgsb1.h
 *
 *  Modified by Vasyl Zhabotynsky on 04/05/2011.
 *
 */
#include <R.h>
#include "Defn1.h"

void lbfgsb1(int n, int m, double *x, double *l, double *u, int *nbd,
	    double *Fmin, optimfn1 fminfn1, optimgr1 fmingr1, int *fail,
	    void *ex, double factr, double pgtol,
	    int *fncount, int *grcount, int maxit, char *msg,
	    int trace, int nREPORT, double *wa, int *iwa, double *g);
