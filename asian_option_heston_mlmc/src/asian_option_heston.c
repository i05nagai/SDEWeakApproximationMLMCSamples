/*
 * $Source: /Users/mariko/lib/source/sde_wa/sde_wa_manual/sde_wa_manual_progra
 * $Revision: 1.6 $
 * $Author: mariko $
 * $Date: 2011/11/29 06:25:53 $
 */
#include <math.h>
#include <sde_wa/sde_wa.h>
#include "asian_option_heston.h"

int ah_V_0(const double y[], double dy[], void *params){
  struct AH_params *pparams;
  pparams=params;
  dy[0]=pparams->mu*y[0];
  dy[1]=pparams->alpha*(pparams->theta-y[1]);
  dy[2]=y[0];
  return SDE_WA_SUCCESS;
}

int ah_str_V_0(const double y[], double dy[], void *params){
  struct AH_params *pparams;
  pparams=params;
  dy[0]=y[0]*(pparams->mu-y[1]/2.0-pparams->rho*pparams->beta/4.0);
  dy[1]=pparams->alpha*(pparams->theta-y[1])-pparams->beta*pparams->beta/4.0;
  dy[2]=y[0];
  return SDE_WA_SUCCESS;
}

int ah_V_1(const double y[], double dy[], void *params){
  struct AH_params *pparams;
  pparams=params;
  dy[0]=y[0]*sqrt(y[1]);
  dy[1]=(pparams->rho)*(pparams->beta)*sqrt(y[1]);
  dy[2]=0.0;
  return SDE_WA_SUCCESS;
}

int ah_V_2(const double y[], double dy[], void *params){
  struct AH_params *pparams;
  pparams=params;
  dy[0]=0.0;
  dy[1]=(pparams->beta)*sqrt((1.0-(pparams->rho)*(pparams->rho))*y[1]);
  dy[2]=0.0;
  return SDE_WA_SUCCESS;
}

int diff_ah_V_1(const double y[], double *dVdy[], void *params){
  struct AH_params *pparams;
  pparams=params;
  dVdy[0][0]=sqrt(y[1]);
  dVdy[0][1]=(y[0]/sqrt(y[1]))/2.0;
  dVdy[0][2]=0.0;
  dVdy[1][0]=0.0;
  dVdy[1][1]=pparams->rho*pparams->beta/(sqrt(y[1])*2.0);
  dVdy[1][2]=0.0;
  dVdy[2][0]=0.0;
  dVdy[2][1]=0.0;
  dVdy[2][2]=0.0;
  return SDE_WA_SUCCESS;
}

int diff_ah_V_2(const double y[], double *dVdy[], void *params){
  struct AH_params *pparams;
  pparams=params;
  dVdy[0][0]=0.0;
  dVdy[0][1]=0.0;
  dVdy[0][2]=0.0;
  dVdy[1][0]=0.0;
  dVdy[1][1]=pparams->beta*sqrt(1.0-pparams->rho*pparams->rho)
    /(sqrt(y[1])*2.0);
  dVdy[1][2]=0.0;
  dVdy[2][0]=0.0;
  dVdy[2][0]=0.0;
  dVdy[2][0]=0.0;
  return SDE_WA_SUCCESS;
}

int exp_ah_sVy_0(double s, const double y[], double dy[], void *params){
double J;
double A;
  struct AH_params *pparams;
  pparams=params;
  J=pparams->theta-pparams->beta*pparams->beta/(4.0*pparams->alpha);
  A=pparams->mu-pparams->rho*pparams->beta*pparams->beta/4.0-y[1]/2.0;
  dy[0]=y[0]*exp((pparams->mu-pparams->rho*pparams->beta*pparams->beta/4.0-J/2.0)*s
                 +(y[1]-J)*(exp(-pparams->alpha*s)-1.0)/(2.0*pparams->alpha));
  dy[1]=J+(y[1]-J)*exp(-pparams->alpha*s);
  dy[2]=y[2]+y[0]*(exp(A*s)-1.0)/A;
  return SDE_WA_SUCCESS;
}

int exp_ah_sVy_1(double s, const double y[], double dy[], void *params){
  struct AH_params *pparams;
  pparams=params;
  dy[0]=y[0]*exp(s*sqrt(y[1])+pparams->rho*pparams->beta*s*s/4.0);
  dy[1]=(pparams->rho*pparams->beta*s/2.0+sqrt(y[1]))*
    (pparams->rho*pparams->beta*s/2.0+sqrt(y[1]));
  dy[2]=y[2];
  return SDE_WA_SUCCESS;
}

int exp_ah_sVy_2(double s, const double y[], double dy[], void *params){
  struct AH_params *pparams;
  pparams=params;
  dy[0]=y[0];
  dy[1]=(s*pparams->beta*sqrt(1.0-pparams->rho*pparams->rho)/2.0 +sqrt(y[1]))
    *(s*pparams->beta*sqrt(1.0-pparams->rho*pparams->rho)/2.0 +sqrt(y[1]));
  dy[2]=y[2];
  return SDE_WA_SUCCESS;
}

double asian_option_heston_call_payoff(double *x, void *params){
	struct AH_params *pparams;
	pparams=params;
	return (x[2] - pparams->K > 0.0) ? (x[2] - pparams->K): 0.0;
}

