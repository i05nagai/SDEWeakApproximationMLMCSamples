/**
 * @file    asian_option_heston.h
 * @brief   header file of asian_option_heston.c.

 *          header file of asian_option_heston.c. 
 * @date    09/Mar/2015 (Mon) 17:32 
 * @author  i05nagai
 */
#ifndef _SDE_WA_MANUAL_PROGRAM_ABS_SYSTEM_H_
#define _SDE_WA_MANUAL_PROGRAM_ABS_SYSTEM_H_
#ifndef M_PIl
#define M_PIl (3.1415926535897932384626433832795029)
#endif
/* parameters for an asian blackscholes model*/
struct ABS_params{
  double mu;
	double sigma;
	double T;
	double K;
};
int abs_V_0(const double y[], double dy[], void *params);
int abs_str_V_0(const double y[], double dy[], void *params);
int abs_V_1(const double y[], double dy[], void *params);
int diff_abs_V_1(const double y[], double *dVdy[], void *params);
int exp_abs_sVy_0(double s, const double y[], double dy[], void *params);
int exp_abs_sVy_1(double s, const double y[], double dy[], void *params);
/* declaration of f : not to be included in SDE_WA_SYSTEM */
double asian_option_blackscholes_call_payoff(double *x, void *params);
#endif
