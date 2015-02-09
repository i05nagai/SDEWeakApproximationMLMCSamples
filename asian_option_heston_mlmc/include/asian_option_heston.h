/*
 * $Source: /Users/mariko/lib/source/sde_wa/sde_wa_manual/sde_wa_manual_progra
 * $Revision: 1.2 $
 * $Author: mariko $
 * $Date: 2011/11/29 05:46:46 $
*/
#ifndef _SDE_WA_MANUAL_PROGRAM_AH_SYSTEM_H_
#define _SDE_WA_MANUAL_PROGRAM_AH_SYSTEM_H_
#ifndef M_PIl
#define M_PIl (3.1415926535897932384626433832795029)
#endif
/* parameters for an asian heston model*/
struct AH_params{
  double alpha;
  double beta;
  double mu;
  double rho;
  double theta;
	double T;
	double K;
};
int ah_V_0(const double y[], double dy[], void *params);
int ah_str_V_0(const double y[], double dy[], void *params);
int ah_V_1(const double y[], double dy[], void *params);
int ah_V_2(const double y[], double dy[], void *params);
int diff_ah_V_1(const double y[], double *dVdy[], void *params);
int diff_ah_V_2(const double y[], double *dVdy[], void *params);
int exp_ah_sVy_0(double s, const double y[], double dy[], void *params);
int exp_ah_sVy_1(double s, const double y[], double dy[], void *params);
int exp_ah_sVy_2(double s, const double y[], double dy[], void *params);

/* declaration of f : not to be included in SDE_WA_SYSTEM */
double asian_option_heston_call_payoff(double *x, void *params);
#endif
