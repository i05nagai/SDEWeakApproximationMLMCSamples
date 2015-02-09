/**
 * @file    european_option_blackscholes.c
 * @brief   Vector fields and the solutions of BlackScholes Model for european options.

 *          Vector fields and the solutions of BlackScholes Model for european options. 
 * @date    09/Mar/2015 (Mon) 17:32 
 * @author  i05nagai
 */
#include <math.h>
#include <sde_wa/sde_wa.h>
#include "european_option_blackscholes.h"

int ebs_V_0(const double y[], double dy[], void *params){
  struct EBS_params *pparams;
  pparams=params;
  dy[0]=pparams->mu*y[0];
  return SDE_WA_SUCCESS;
}

int ebs_str_V_0(const double y[], double dy[], void *params){
  struct EBS_params *pparams;
  pparams=params;
  dy[0]=y[0]*(pparams->mu - pparams->sigma*pparams->sigma/2.0);
  return SDE_WA_SUCCESS;
}

int ebs_V_1(const double y[], double dy[], void *params){
  struct EBS_params *pparams;
  pparams=params;
  dy[0]=y[0]*pparams->sigma;
  return SDE_WA_SUCCESS;
}

int diff_ebs_V_1(const double y[], double *dVdy[], void *params){
  struct EBS_params *pparams;
  pparams=params;
  dVdy[0][0]=pparams->sigma;
  return SDE_WA_SUCCESS;
}

int exp_ebs_sVy_0(double s, const double y[], double dy[], void *params){
  struct EBS_params *pparams;
  pparams=params;
  dy[0]=y[0]*exp((pparams->mu - pparams->sigma*pparams->sigma/2.0)*s);
  return SDE_WA_SUCCESS;
}

int exp_ebs_sVy_1(double s, const double y[], double dy[], void *params){
  struct EBS_params *pparams;
  pparams=params;
  dy[0]=y[0]*exp(s*pparams->sigma);
  return SDE_WA_SUCCESS;
}

double european_option_blackscholes_call_payoff(double *x, void *params){
	struct EBS_params *pparams;
	pparams=params;
	return (x[0] - pparams->K > 0.0) ? (x[0] - pparams->K): 0.0;
}


