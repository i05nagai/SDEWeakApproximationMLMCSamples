/**
 * @file    asian_option_heston_mlmc_nv.c
 * @brief   level l estimation for the Ninomiya-Victoir.

 *          level l estimation for the Ninomiya-Victoir. 
 * @date    09/Mar/2015 (Mon) 17:32 
 * @author  i05nagai
 */
#include <stdlib.h> /* for malloc */
#include <stdio.h>
#include <math.h>
#include <sde_wa/sde_wa.h>
#include <sde_wa/sde_wa_em.h>
#include "asian_option_heston.h"
#include "asian_option_heston_mlmc.h"
#include <sde_mlmc/mlmc.h>
#include <mt19937/mt19937.h>


/**
 * level l estimaton
 * level l esitmator using Ninomiya Victoir
 * @param mlmc 		a poinster of SDE_MLMC.
 * @param estimator	a pointer of MLMC_L.
 * @param l			level
 * @param M			a base of the number of partitions 
 * @param N			the num of samles for level l estimation.
 * @return None.
 */
int mlmc_level_l_estimation_nv(SDE_MLMC *mlmc, MLMC_L *estimator, int l, int M, int N) {
	SDE_WA_SLTN *sl = mlmc->sl;
	SDE_WA_SYSTEM *sde = sl->sde;
	int i,j,k,m;
	int n = pow(M, l);
	double dt1, dt2;
	double *u_seq;
	double *n_seq;
	double *sp;
	double *tmp_pt;
	double tmp;
	int neg_except_flag=0;
	int neg_event_counter=0;
	double **xf;
	double **xc;

	RV_NV *rv_nv;
	rv_nv=(RV_NV *)malloc(sizeof(RV_NV));
	rv_nv->rv_nv_n=(double *)malloc(sizeof(double)*sde->dim_BM);

	xf=(double **)malloc(sizeof(double *)*2);
	for (i=0; i<2; i++) xf[i]=(double *)malloc(sizeof(double)*sde->dim_y);

	xc=(double **)malloc(sizeof(double *)*2);
	for (i=0; i<2; i++) xc[i]=(double *)malloc(sizeof(double)*sde->dim_y);

	int n_dim = n*sde->dim_BM;
	int rv_dim = n_dim + n_dim%2 + n;
	u_seq=(double *)malloc(sizeof(double)*(rv_dim));
	n_seq=(double *)malloc(sizeof(double)*(n_dim));
	sp=(double *)malloc(sizeof(double)*sde->dim_BM);

	if (l == 0) {
		for (j=0; j<sde->dim_y; j++) {
			xf[0][j] = mlmc->init_y[j];
		}
		for (neg_except_flag=0, i=0; i < N; i++){
			for (j=0; j<rv_dim; j++) {
				u_seq[j] = genrand_real2();
			}
			trans_rand_unif_to_normal(u_seq, rv_nv->rv_nv_n, sde->dim_BM);

			rv_nv->rv_nv_b=(int)(u_seq[n_dim + n_dim%2]+0.5);
			next_SDE_WA(sl, 1.0, xf[0], xf[1], rv_nv);
			if (xf[1][1] < 0.0 || isnan(xf[1][1]) || isnan(xf[1][0])) {
				i-=1; neg_except_flag=0; neg_event_counter++;
				continue;
			}

			estimator->sumY += mlmc->payoff_function(xf[1], sde->params);
			estimator->sumYY += mlmc->payoff_function(xf[1], sde->params) * mlmc->payoff_function(xf[1], sde->params);
		} /* for i */
	} else {
		for (neg_except_flag=0, i=0; i < N; i++){
			for (j=0; j<rv_dim; j++) {
				u_seq[j] = genrand_real2();
			}
			trans_rand_unif_to_normal(u_seq, n_seq, n_dim);

			dt1=mlmc->T/(double)n;
			dt2=mlmc->T*M/(double)n;
			for (j=0; j<sde->dim_y; j++) {
				xf[0][j] = mlmc->init_y[j];
				xc[0][j] = mlmc->init_y[j];
			}
			for (k=0; k<n/M; k++){
				for (m=0; m<sde->dim_BM; m++) {
					sp[m] = 0.0;
				}
				for (j=0; j<M; j++){
					sp[m] += rv_nv->rv_nv_n[m];
					for (m=0; m<sde->dim_BM; m++) {
						rv_nv->rv_nv_n[m] = n_seq[k*M*sde->dim_BM + m + sde->dim_BM*j];
						sp[m] += rv_nv->rv_nv_n[m];
					}
					rv_nv->rv_nv_b = (int)(u_seq[n_dim+n_dim%2 + k*M + j]+0.5);
					next_SDE_WA(sl, dt1, xf[0], xf[1], rv_nv);
					if (xf[1][1] < 0.0 || isnan(xf[1][1]) || isnan(xf[1][0])) {
						neg_except_flag=1; break;
					}
					tmp_pt=xf[0];
					xf[0]=xf[1];
					xf[1]=tmp_pt;
				}
				for (m=0; m<sde->dim_BM; m++) {
					rv_nv->rv_nv_n[m] = sp[m]/sqrt(M);
				}
				next_SDE_WA(sl, dt2, xc[0], xc[1], rv_nv);
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
	free(n_seq);
	free(u_seq);
	free(sp);
	free(rv_nv->rv_nv_n);
	free(rv_nv);
	for (i=0; i<2; i++) {
		free(xf[i]);
		free(xc[i]);
	}
	free(xf);
	free(xc);
	return 0;
}
