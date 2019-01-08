/*
 *  mle_bb.h
 *
 *  Created by Wei Sun on 12/02/2016.
 *
 */

#include <stdio.h>
#include <stddef.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
#include "utility.h"
#include "lbfgsb1.h"

double negLogL (int n, double* para, void* ex);

void negGradLogL (int n, double* para, double* gr, void* ex);

void mle_bb (int* dims, double* nAR, double* nTotalR, double* wsR,
             double* pi, double* rho, double* min_pi, double* max_pi,
             double* min_rho, double* max_rho, double* logLik);
