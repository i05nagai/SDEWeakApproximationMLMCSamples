/*
 * $Source: /Users/mariko/lib/source/sde_wa/sde_wa_manual/sde_wa_manual_progra
 * $Revision: 1.6 $
 * $Author: mariko $
 * $Date: 2011/11/29 06:25:53 $
 */
#include <math.h>
#include <sde_wa/sde_wa.h>
#include "european_option_heston.h"

int eh_V_0(const double y[], double dy[], void *params){
  struct EH_params *pparams;
  pparams=params;
  dy[0]=pparams->mu*y[0];
  dy[1]=pparams->alpha*(pparams->theta-y[1]);
  return SDE_WA_SUCCESS;
}

int eh_str_V_0(const double y[], double dy[], void *params){
  struct EH_params *pparams;
  pparams=params;
  dy[0]=y[0]*(pparams->mu-y[1]/2.0-pparams->rho*pparams->beta/4.0);
  dy[1]=pparams->alpha*(pparams->theta-y[1])-pparams->beta*pparams->beta/4.0;
  return SDE_WA_SUCCESS;
}

int eh_V_1(const double y[], double dy[], void *params){
  struct EH_params *pparams;
  pparams=params;
  dy[0]=y[0]*sqrt(y[1]);
  dy[1]=(pparams->rho)*(pparams->beta)*sqrt(y[1]);
  return SDE_WA_SUCCESS;
}

int eh_V_2(const double y[], double dy[], void *params){
  struct EH_params *pparams;
  pparams=params;
  dy[0]=0.0;
  dy[1]=(pparams->beta)*sqrt((1.0-(pparams->rho)*(pparams->rho))*y[1]);
  return SDE_WA_SUCCESS;
}

int diff_eh_V_1(const double y[], double *dVdy[], void *params){
  struct EH_params *pparams;
  pparams=params;
  dVdy[0][0]=sqrt(y[1]);
  dVdy[0][1]=(y[0]/sqrt(y[1]))/2.0;
  dVdy[1][0]=0.0;
  dVdy[1][1]=pparams->rho*pparams->beta/(sqrt(y[1])*2.0);
  return SDE_WA_SUCCESS;
}

int diff_eh_V_2(const double y[], double *dVdy[], void *params){
  struct EH_params *pparams;
  pparams=params;
  dVdy[0][0]=0.0;
  dVdy[0][1]=0.0;
  dVdy[1][0]=0.0;
  dVdy[1][1]=pparams->beta*sqrt(1.0-pparams->rho*pparams->rho)
    /(sqrt(y[1])*2.0);
  return SDE_WA_SUCCESS;
}

int exp_eh_sVy_0(double s, const double y[], double dy[], void *params){
double J;
double A;
  struct EH_params *pparams;
  pparams=params;
  J=pparams->theta-pparams->beta*pparams->beta/(4.0*pparams->alpha);
  A=pparams->mu-pparams->rho*pparams->beta*pparams->beta/4.0-y[1]/2.0;
  dy[0]=y[0]*exp((pparams->mu-pparams->rho*pparams->beta/4.0-J/2.0)*s
                 +(y[1]-J)*(exp(-pparams->alpha*s)-1.0)/(2.0*pparams->alpha));
  dy[1]=J+(y[1]-J)*exp(-pparams->alpha*s);
  return SDE_WA_SUCCESS;
}

int exp_eh_sVy_00(double s, const double y[], double dy[], void *params){
  struct EH_params *pparams;
  pparams=params;
  dy[0]=y[0]*exp(pparams->mu*s);
  dy[1]=y[1]*exp(-pparams->alpha*s);
  return SDE_WA_SUCCESS;
}

int exp_eh_sVy_01(double s, const double y[], double dy[], void *params){
  struct EH_params *pparams;
  pparams=params;
  double J=pparams->alpha*pparams->theta - pparams->beta*pparams->beta/4.0;
  dy[0]=y[0]*exp(-(J*s*s/4.0 + (y[1]/2.0 + pparams->rho*pparams->beta/4.0)*s));
  dy[1]=y[1] + J*s;
  return SDE_WA_SUCCESS;
}

int exp_eh_sVy_1(double s, const double y[], double dy[], void *params){
  struct EH_params *pparams;
  pparams=params;
  dy[0]=y[0]*exp(s*sqrt(y[1])+pparams->rho*pparams->beta*s*s/4.0);
  dy[1]=(pparams->rho*pparams->beta*s/2.0+sqrt(y[1]))*
    (pparams->rho*pparams->beta*s/2.0+sqrt(y[1]));
  return SDE_WA_SUCCESS;
}

int exp_eh_sVy_2(double s, const double y[], double dy[], void *params){
  struct EH_params *pparams;
  pparams=params;
  dy[0]=y[0];
  dy[1]=(s*pparams->beta*sqrt(1.0-pparams->rho*pparams->rho)/2.0 +sqrt(y[1]))
    *(s*pparams->beta*sqrt(1.0-pparams->rho*pparams->rho)/2.0 +sqrt(y[1]));
  return SDE_WA_SUCCESS;
}

double european_option_heston_call_payoff(double *x, void *params){
	struct EH_params *pparams;
	pparams=params;
	return (x[0] - pparams->K > 0.0) ? (x[0] - pparams->K): 0.0;
}

