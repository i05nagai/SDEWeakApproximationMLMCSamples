/**
 * @file    european_option_blackscholes_mlmc_em.c
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
#include "european_option_blackscholes.h"
#include "european_option_blackscholes_mlmc.h"
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
	double **xf;
	double **xc;

	u_seq=(double *)malloc(sizeof(double)*n*sde->dim_BM);
	n_seq=(double *)malloc(sizeof(double)*n*sde->dim_BM);

	xf=(double **)malloc(sizeof(double *)*2);
	for (i=0; i<2; i++) xf[i]=(double *)malloc(sizeof(double)*sde->dim_y);

	xc=(double **)malloc(sizeof(double *)*2);
	for (i=0; i<2; i++) xc[i]=(double *)malloc(sizeof(double)*sde->dim_y);

	if (l == 0) {
		sp=(double *)malloc(sizeof(double)*sde->dim_BM*2);
		sp2=(double *)malloc(sizeof(double)*sde->dim_BM*2);
		for (j=0; j<sde->dim_y; j++) {
			xf[0][j] = mlmc->init_y[j];
		}
		for (i=0; i < N; i++){
			for (j=0; j<(sde->dim_BM + sde->dim_BM%2); j++) {
				u_seq[j] = genrand_real2();
			}
			trans_rand_unif_to_normal(u_seq, n_seq, sde->dim_BM*2);

			next_SDE_WA(sl, mlmc->T, xf[0], xf[1], n_seq);

			estimator->sumY += mlmc->payoff_function(xf[1], sde->params);
			estimator->sumYY += mlmc->payoff_function(xf[1], sde->params) * mlmc->payoff_function(xf[1], sde->params);
		} /* for i */
	} else {
		sp=(double *)malloc(sizeof(double)*sde->dim_BM);
		sp2=(double *)malloc(sizeof(double)*sde->dim_BM);
		dt1=mlmc->T/(double)n;
		dt2=mlmc->T*M/(double)n;
		for (i=0; i < N; i++){
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
						sp[m] = n_seq[k*M*sde->dim_BM + m + sde->dim_BM*j];
						sp2[m] += sp[m];
					}
					next_SDE_WA(sl, dt1, xf[0], xf[1], sp);
					tmp_pt=xf[0];
					xf[0]=xf[1];
					xf[1]=tmp_pt;
				}
				for (m=0; m<sde->dim_BM; m++) {
					sp2[m] = sp2[m]/sqrt(M);
				}
				next_SDE_WA(sl, dt2, xc[0], xc[1], sp2);
				tmp_pt=xc[0];
				xc[0]=xc[1];
				xc[1]=tmp_pt;
			}
			tmp = mlmc->payoff_function(xf[0], sde->params) - mlmc->payoff_function(xc[0], sde->params);
			estimator->sumY += tmp; 
			estimator->sumYY += tmp*tmp;
		} /* for i */
	}
	return 0;
}

