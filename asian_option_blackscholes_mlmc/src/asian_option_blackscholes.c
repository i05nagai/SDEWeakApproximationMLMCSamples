/**
 * @file    asian_option_blackscholes.c
 * @brief   Vector fields and the solutions of BlackScholes Model for asian options.

 *          Vector fields and the solutions of BlackScholes Model for asian options. 
 * @date    09/Mar/2015 (Mon) 17:32 
 * @author  i05nagai
 */
#include <math.h>
#include <sde_wa/sde_wa.h>
#include "asian_option_blackscholes.h"

int abs_V_0(const double y[], double dy[], void *params){
  struct ABS_params *pparams;
  pparams=params;
  dy[0]=pparams->mu*y[0];
  dy[1]=y[0];
  return SDE_WA_SUCCESS;
}

int abs_str_V_0(const double y[], double dy[], void *params){
  struct ABS_params *pparams;
  pparams=params;
  dy[0]=y[0]*(pparams->mu-pparams->sigma*pparams->sigma/2.0);
  dy[1]=y[0];
  return SDE_WA_SUCCESS;
}

int abs_V_1(const double y[], double dy[], void *params){
  struct ABS_params *pparams;
  pparams=params;
  dy[0]=y[0]*pparams->sigma;
  dy[1]=0.0;
  return SDE_WA_SUCCESS;
}

int diff_abs_V_1(const double y[], double *dVdy[], void *params){
  struct ABS_params *pparams;
  pparams=params;
  dVdy[0][0]=pparams->sigma;
  dVdy[0][1]=0.0;
  dVdy[1][0]=0.0;
  dVdy[1][1]=0.0;
  return SDE_WA_SUCCESS;
}

int exp_abs_sVy_0(double s, const double y[], double dy[], void *params){
  struct ABS_params *pparams;
  pparams=params;
double J = (pparams->mu - pparams->sigma*pparams->sigma/2.0);
  dy[0]=y[0]*exp(s*J);
  dy[1]=y[1] + y[0]*(exp(s*J)-1.0)/J;
  return SDE_WA_SUCCESS;
}

int exp_abs_sVy_1(double s, const double y[], double dy[], void *params){
  struct ABS_params *pparams;
  pparams=params;
  dy[0]=y[0]*exp(s*pparams->sigma);
  dy[1]=y[1];
  return SDE_WA_SUCCESS;
}

double asian_option_blackscholes_call_payoff(double *x, void *params){
	struct ABS_params *pparams;
	pparams=params;
	return (x[1] - pparams->K > 0.0) ? (x[1] - pparams->K): 0.0;
}


