/**
 * @file    european_option_heston_mlmc_em.c
 * @brief   level l estimation for the Euler Maruyama Scheme.

 *          level l estimation for the Euler Maruyama Scheme. 
 * @date    09/Mar/2015 (Mon) 17:32 
 * @author  i05nagai
 */
#include <stdlib.h> /* for malloc */
#include <stdio.h>
#include <math.h>
#include <sde_wa/sde_wa.h>
#include <sde_wa/sde_wa_em.h>
#include "european_option_heston.h"
#include "european_option_heston_mlmc.h"
#include <sde_mlmc/mlmc.h>
#include <mt19937/mt19937.h>

/**
 * level l estimaton
 * level l esitmator using Euler Maruyama 
 * @param mlmc 		a poinster of SDE_MLMC.
 * @param estimator	a pointer of MLMC_L.
 * @param l			level
 * @param M			a base of the number of partitions 
 * @param N			the num of samles for level l estimation.
 * @return None.
 */
int mlmc_level_l_estimation_em(SDE_MLMC *mlmc, MLMC_L *estimator, int l, int M, int N) {
	SDE_WA_SLTN *sl = mlmc->sl;
	SDE_WA_SYSTEM *sde = sl->sde;
	int i,j,k,m;
	int n = pow(M, l);
	double dt1, dt2;
	double *u_seq;
	double *n_seq;
	double *sp, *sp2;
	double *tmp_pt;
	double tmp;
	int neg_event_counter=0;
	int neg_except_flag=0;
	//struct EH_params *params = sde->params;
	double **xf;
	double **xc;

	u_seq=(double *)malloc(sizeof(double)*n*sde->dim_BM);
	n_seq=(double *)malloc(sizeof(double)*n*sde->dim_BM);
	sp=(double *)malloc(sizeof(double)*sde->dim_BM);
	sp2=(double *)malloc(sizeof(double)*sde->dim_BM);

	xf=(double **)malloc(sizeof(double *)*2);
	for (i=0; i<2; i++) xf[i]=(double *)malloc(sizeof(double)*sde->dim_y);

	xc=(double **)malloc(sizeof(double *)*2);
	for (i=0; i<2; i++) xc[i]=(double *)malloc(sizeof(double)*sde->dim_y);

	if (l == 0) {
		for (j=0; j<sde->dim_y; j++) {
			xf[0][j] = mlmc->init_y[j];
		}
		for (neg_except_flag=0, i=0; i < N; i++){
			for (j=0; j<(sde->dim_BM + sde->dim_BM%2); j++) {
				u_seq[j] = genrand_real2();
			}
			trans_rand_unif_to_normal(u_seq, n_seq, sde->dim_BM);

			next_SDE_WA(sl, mlmc->T, xf[0], xf[1], n_seq);
			//xf[1][1] = xf[0][1] + (1.0 - exp(-params->alpha*dt1))*((params->theta - xf[0][1])) + params->beta*0.25*sqrt(xf[0][1]*dt1)*(params->rho*sp[0] + sqrt(1-params->rho*params->rho)*sp[1]);
			if (xf[1][1] < 0.0 || isnan(xf[1][1]) || isnan(xf[1][0])) {
				i-=1; neg_event_counter++;
				continue;
			}
			tmp =mlmc->payoff_function(xf[1], sde->params);
			estimator->sumY += tmp;
			estimator->sumYY += tmp*tmp;
		} /* for i */
	} else {
		dt1=mlmc->T/(double)n;
		dt2=mlmc->T*M/(double)n;
		for (neg_except_flag=0, i=0; i < N; i++){
			for (j=0; j<n*sde->dim_BM+n*sde->dim_BM%2; j++) {
				u_seq[j] = genrand_real2();
			}
			trans_rand_unif_to_normal(u_seq, n_seq, n*sde->dim_BM);

			for (j=0; j<sde->dim_y; j++) {
				xf[0][j] = mlmc->init_y[j];
				xc[0][j] = mlmc->init_y[j];
			}
			for (k=0; k<n/M; k++){
				for (m=0; m<sde->dim_BM; m++) {
					sp2[m] = 0.0;
				}
				for (j=0; j<M; j++){
					for (m=0; m<sde->dim_BM; m++) {
						sp[m] = n_seq[k*M + m + sde->dim_BM*j];
						sp2[m] += sp[m];
					}
					next_SDE_WA(sl, dt1, xf[0], xf[1], sp);
					//xf[1][1] = params->theta - exp(-params->alpha*dt1)*((params->theta - xf[0][1])) + params->beta*0.25*sqrt(xf[0][1]*dt1)*(params->rho*sp[0] + sqrt(1-params->rho*params->rho)*sp[1]);

					if (xf[1][1] < 0.0 || isnan(xf[1][1]) || isnan(xf[1][0])) {
						neg_except_flag=1; break;
					}
					tmp_pt=xf[0];
					xf[0]=xf[1];
					xf[1]=tmp_pt;
				}
				for (m=0; m<sde->dim_BM; m++) {
					sp2[m] = sp2[m]/sqrt(M);
				}
				next_SDE_WA(sl, dt2, xc[0], xc[1], sp2);
				//xc[1][1] = params->theta - exp(-params->alpha*dt2)*((params->theta - xc[0][1])) + params->beta*0.25*sqrt(xc[0][1]*dt2)*(params->rho*sp2[0] + sqrt(1-params->rho*params->rho)*sp2[1]);
				if (xc[1][1] < 0.0 || isnan(xc[1][1]) || isnan(xc[1][0])) {
					neg_except_flag=1; break;
				}
				tmp_pt=xc[0];
				xc[0]=xc[1];
				xc[1]=tmp_pt;
			}
			if (neg_except_flag) {
				i-=1; neg_except_flag=0; neg_event_counter++;
				continue;
			}
			tmp = mlmc->payoff_function(xf[0], sde->params) - mlmc->payoff_function(xc[0], sde->params);
			estimator->sumY += tmp; 
			estimator->sumYY += tmp*tmp;
		} /* for i */
	}
	printf("neg_event_counter:%d\n", neg_event_counter);
	return 0;
}

