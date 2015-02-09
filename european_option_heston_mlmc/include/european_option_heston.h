/*
 * $Source: /Users/mariko/lib/source/sde_wa/sde_wa_manual/sde_wa_manual_progra
 * $Revision: 1.2 $
 * $Author: mariko $
 * $Date: 2011/11/29 05:46:46 $
*/
#ifndef _SDE_WA_MANUAL_PROGRAM_EH_SYSTEM_H_
#define _SDE_WA_MANUAL_PROGRAM_EH_SYSTEM_H_
#ifndef M_PIl
#define M_PIl (3.1415926535897932384626433832795029)
#endif
struct EH_params{
  double alpha;
  double beta;
  double mu;
  double rho;
  double theta;
	double T;
	double K;
};
int eh_V_0(const double y[], double dy[], void *params);
int eh_str_V_0(const double y[], double dy[], void *params);
int eh_V_1(const double y[], double dy[], void *params);
int eh_V_2(const double y[], double dy[], void *params);
int diff_eh_V_1(const double y[], double *dVdy[], void *params);
int diff_eh_V_2(const double y[], double *dVdy[], void *params);
int exp_eh_sVy_0(double s, const double y[], double dy[], void *params);
int exp_eh_sVy_00(double s, const double y[], double dy[], void *params);
int exp_eh_sVy_01(double s, const double y[], double dy[], void *params);
int exp_eh_sVy_1(double s, const double y[], double dy[], void *params);
int exp_eh_sVy_2(double s, const double y[], double dy[], void *params);
/* declaration of f : not to be included in SDE_WA_SYSTEM */
double european_option_heston_call_payoff(double *x, void *params);
#endif
