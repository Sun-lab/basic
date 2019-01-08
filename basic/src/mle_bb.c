/*
 *  mle_bb.c
 *
 *  Created by Wei Sun on 12/02/2016.
 *
 */

#include "mle_bb.h"

/**********************************************************************
 *
 * negative log likelihood and gradient function
 *
 * The first parameter (n) is the number of paramters, n=2
 *
 * Sample size N is included in the parameter "ex".
 *
 **********************************************************************/

double negLogL (int n, double* para, void* ex){
  int i, k, N;
  double sumL, logLi, ni, ni0;
  double pi, rho, theta;
  double *exPara, *nA, *nTotal, *ws, *lchooseSaved;
  
  exPara = (double *) ex;
  N      = ceil(exPara[0]-0.5);

  nA     = exPara + 1;
  nTotal = nA + N;
  ws     = nTotal + N;
  lchooseSaved = ws + N;

  pi    = para[0];
  rho   = para[1];
  theta = rho/(1 - rho);

  sumL = 0.0;
  
  for(i=0; i<N; i++){
    ni   = nTotal[i];
    ni0  = nA[i];
    
    logLi = lchooseSaved[i];
    
    if(ni0 > 0){
      for(k=0; k<ni0; k++) logLi += log(pi + k*theta);
    }
    
    if(ni0 < ni){
      for(k=0; k<ni-ni0; k++) logLi += log(1.0 - pi + k*theta);
    }
    
    for(k=0; k<ni; k++){
      logLi -= log(1.0 + k*theta);
    }
    
    sumL = sumL + logLi*ws[i];
    
  }
  
  return(-sumL);
}

/**********************************************************************
 *
 * negGradLogH1: gradient negLogH1
 *
 **********************************************************************/

void negGradLogL (int n, double* para, double* gr, void* ex)
{
  int i, k, N;
  double gradPi, gradTh, gradPi_i, gradTh_i, pi, rho, theta, tmp;
  double *exPara, *nA, *nTotal, *ws, ni, ni0;
  
  exPara = (double *) ex;
  N      = ceil(exPara[0]-0.5);

  nA     = exPara + 1;
  nTotal = nA + N;
  ws     = nTotal + N;
  
  pi     = para[0];
  rho    = para[1];
  theta  = rho/(1 - rho);

  gradPi = 0.0;
  gradTh = 0.0;
  
  for(i=0; i<N; i++){
    
    ni   = nTotal[i];
    ni0  = nA[i];
    
    gradPi_i = 0.0;
    gradTh_i = 0.0;

    if(ni0 > 0){
      for(k=0; k<ni0; k++){
        tmp = 1.0/(pi + k*theta);
        
        gradPi_i += tmp;
        gradTh_i += k*tmp;
      }
    }
    
    if(ni0 < ni){
      for(k=0; k<ni-ni0; k++){
        tmp = 1.0/(1.0 - pi + k*theta);
        
        gradPi_i -= tmp;
        gradTh_i += k*tmp;
      }
    }
    
    for(k=0; k<ni; k++){
      gradTh_i -=  k/(1 + k*theta);
    }

    gradPi = gradPi + gradPi_i*ws[i];
    gradTh = gradTh + gradTh_i*ws[i];

  }
  
  gr[0] = -gradPi;
  gr[1] = -gradTh/(1-rho)/(1-rho);
}

/**********************************************************************
 *
 * mle_bb
 *
 * identify (weighted) MLE of beta-binomial distribution
 *
 **********************************************************************/


void mle_bb (int* dims, double* nAR, double* nTotalR, double* wsR,
             double* pi, double* rho, double* min_pi, double* max_pi,
             double* min_rho, double* max_rho, double* logLik){
  int k;
  double *exPara, *nA, *nTotal, *ws, *lchooseSaved;
  int N = dims[0]; //sample size
  
  /** 
   * we have to combine nA, nTotal, and zeta into one vector for
   * its usage in sovling for MLE
   */
  
  exPara = (double *) R_alloc(4*N+1, sizeof(double));
  nA     = exPara + 1;
  nTotal = nA + N;
  ws     = nTotal + N;
  lchooseSaved = ws + N;

  /**
   * exPara[0]  is sample size, fixed throught the computation 
   */
  exPara[0] = (double) N; 

  
  for (k=0; k<N; k++) {
    nTotal[k] = nTotalR[k];
    nA[k]     = nAR[k];
    ws[k]     = wsR[k];
    // these lchoose constant will be calculated many times
    // during the optimzation, so we pre-compute them here
    lchooseSaved[k] = lchoose(nTotal[k], nA[k]);
  }
  
  
  /* **********************************************************
   * parameters for function lbfgsb, which will be used to
     obtain MLE for H1: with allelic imbalance
   
   void lbfgsb(int n, int lmm, double *x, double *lower,
          double *upper, int *nbd, double *Fmin, optimfn fn,
          optimgr gr, int *fail, void *ex, double factr,
          double pgtol, int *fncount, int *grcount,
          int maxit, char *msg, int trace, int nREPORT);
   
   n:       the number of parameters
   
   lmm:     is an integer giving the number of BFGS updates 
            retained in the "L-BFGS-B" method, It defaults to 5.
   
   x:       starting parameters on entry and the final parameters on exit
   
   lower:   lower bounds
   
   upper:   upper bounds
   
   nbd:     specifies which bounds are to be used. 
            nbd(i)=0 if x(i) is unbounded,
            1 if x(i) has only a lower bound,
            2 if x(i) has both lower and upper bounds, and
            3 if x(i) has only an upper bound.
            On exit nbd is unchanged.
   
   Fmin:    final value of the function
   
   fn:      the function to be minimized
   
   gr:      the gradient function
   
   fail:    integer code, 0 for success, 51 for warning and 52 for error
   
   ex:      extra parameters for the function to be minimized
   
   factr:   controls the convergence of the "L-BFGS-B" method. 
            Convergence occurs when the reduction in the objective is 
            within this factor of the machine tolerance. Default is 1e7, 
            that is a tolerance of about 1e-8.
   
   pgtol:   helps control the convergence of the "L-BFGS-B" method. 
            It is a tolerance on the projected gradient in the current 
            search direction. This defaults to zero, when the check 
            is suppressed.
   
   fncount: the number of calls to fn 
   
   grcount: the number of calls to gr
   
   maxit:   maximum of iterations
   
   msg:     A character string giving any additional information 
            returned by the optimizer, or NULL
   
   trace:   Non-negative integer. If positive, tracing information 
            on the progress of the optimization is produced. 
            Higher values may produce more tracing information: 
            for method "L-BFGS-B" there are six levels of tracing. 
   
   nREPORT: The frequency of reports for the "BFGS", "L-BFGS-B" 
            and "SANN" methods if control$trace is positive. 
            Defaults to every 10 iterations for "BFGS" and "L-BFGS-B"
   
   * **********************************************************/
  
  int npara, lmm, fail, fncount, grcount, maxit, nREPORT;
  int nbd[2];

  npara   = 2;
  lmm     = 5;
  fail    = 0;
  fncount = 0;
  grcount = 0;
  maxit   = 100;
  nREPORT = 5;
  
  nbd[0]  = 2;
  nbd[1]  = 2;
  
  //technical parameters below:
  double *wa, *g1;
  int *iwa;
  
  //consider replacing with simple Calloc
  wa  = (double *) S_alloc(2*lmm*npara+4*npara+11*lmm*lmm+8*lmm,sizeof(double));
  iwa = (int*) R_alloc(3*npara,sizeof(int));
  g1  = (double *)R_alloc(npara, sizeof(double));

  double initPara[2];
  double lower[2];
  double upper[2];
  double Fmin, factr, pgtol;
  
  /* initPara = c(pi, rho) */
  
  initPara[0] = *pi;
  initPara[1] = *rho;

  lower[0] = *min_pi;
  upper[0] = *max_pi;
  
  lower[1] = *min_rho;
  upper[1] = *max_rho;

  factr = 1e7;
  pgtol = 0.0;
  
  char msg[1023];
  
  lbfgsb1(npara, lmm, initPara, lower, upper, nbd, &Fmin,
         negLogL, negGradLogL, &fail, (void*)exPara, factr, pgtol,
         &fncount, &grcount, maxit, msg, 0, nREPORT, wa, iwa, g1);

  if (fail) {
    // Rprintf("not a clearn success of L-BFGS-B optimization: %d\n", fail);
    // Rprintf("%s \n", msg);
    
  }else{
    *pi  = initPara[0];
    *rho = initPara[1];
    *logLik = -Fmin;
  }

}
