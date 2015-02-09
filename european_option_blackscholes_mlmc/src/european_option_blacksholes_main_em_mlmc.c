/*
 * $Source: /Users/mariko/lib/source/sde_wa/sde_wa_manual/sde_wa_manual_progra
 * $Revision: 1.7 $
 * $Author: mariko $
 * $Date: 2011/11/29 10:48:07 $
 */
#include <stdlib.h> /* for malloc */
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_qrng.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <sde_wa.h>
#include "european_option_blacksholes.h"
#include "mlmc.h"

#define L_MAX 20

int main(void){
	enum ALG alg = E_M;
	int m_moment = 7; /* Romberg */
	int M = 4;
	int n;
	double x0=1.0;
	double K = 1.00;
	struct EBS_params ebs_params;
	double mu=0.05;
	double sigma=0.2;
	double T=1.0;
	ebs_params.mu=mu;
	ebs_params.sigma=sigma;
	ebs_params.T=T;
	ebs_params.K=K;
	MLMC_L **mlmc_est;
	{
		SDE_WA_SYSTEM *sde;
		SDE_WA_SLTN *sl;
		SDE_WA_PRICER *pricer;
		double sum;
		int dNl;
		double hl;
		double eps=0.0001;
		double epss[] = {0.001, 0.0005, 0.0002, 0.0001, 0.00005};
		int l;
		int L;
		int i;
		double con;
		int converged_flag=0;
		int extraporation_flag=0;
		sde=alloc_SDE_WA_SYSTEM(1, 1, &ebs_params);
		sde->sde_type=STR;

		sde->V[0]=ebs_str_V_0;
		sde->V[1]=ebs_V_1;

		sde->drift_corrector[0]=NULL;
		sde->drift_corrector[1]=diff_ebs_V_1; /* for EM */

		sde->exp_sV[0]=exp_ebs_sVy_0; /* for NV */
		sde->exp_sV[1]=exp_ebs_sVy_1; /* for NV */

		sl=alloc_SDE_WA_SLTN(alg, m_moment, sde);

		pricer=alloc_SDE_WA_PRICER(sl);
		pricer->payoff_function=european_option_blacksholes_call_payoff;
		pricer->init_y[0] = x0;
		pricer->T = 1.0;

		//gsl_qrng *q;
		gsl_rng_type *rng_type = (gsl_rng_type *)gsl_rng_default;
		gsl_rng *rng = gsl_rng_alloc(rng_type);

		if (pricer->gen_unif_rand == NULL) {
			/* sobol seq. => seq.dim.<=40*/
		//	q=gsl_qrng_alloc(gsl_qrng_sobol, n*sde->dim_BM);
		//	pricer->gen_unif_rand = q;
		//	pricer->next_rand = gen_unif_rand_number_gsl_quasi;

			pricer->gen_unif_rand = rng;
			pricer->next_rand = gen_unif_rand_number_gsl_pseudo;
		}

		mlmc_est = (MLMC_L **)malloc(sizeof(MLMC_L *)*L_MAX);
		for (l=0; l<L_MAX; l++) {
			mlmc_est[l] = alloc_level_l_estimator();
		}

		for (l=0; l<=-4; l++) {
			L=l;
			mlmc_est[L]->N = 4000000;
			mlmc_level_l_estimation(pricer, mlmc_est[L], L, M, mlmc_est[L]->N);
			mlmc_est[L]->V = mlmc_est[L]->sumV/mlmc_est[L]->N - (mlmc_est[L]->sumY/mlmc_est[L]->N)*(mlmc_est[L]->sumY/mlmc_est[L]->N);
			printf("L:%d M:%d\n", L, M);
			printf("N:%d, Y:%12e, V:%12.12lf, sumV:%lf\n", mlmc_est[L]->N, mlmc_est[L]->sumY/mlmc_est[L]->N, mlmc_est[L]->V, mlmc_est[L]->sumV);
		}
		//exit(-1);
		for (i=0; i<5; i++) {
			L=-1;
			eps = epss[i];
			for (l=0; l<L_MAX; l++) {
				init_level_l_estimator(mlmc_est[l]);
			}
			converged_flag=0;
			printf("eps:%lf \n", eps);
			while (!converged_flag) {
				L++;
				n = (int)pow(M,L);

				mlmc_level_l_estimation(pricer, mlmc_est[L], L, M, mlmc_est[L]->N);
				printf("N:%d, Y:%lf, sumV:%lf\n", mlmc_est[L]->N, mlmc_est[L]->sumY/mlmc_est[L]->N, mlmc_est[L]->sumV);
				//estimate Vl
				sum = 0.0;
				for (l=0; l<=L; l++) {
					mlmc_est[l]->V = mlmc_est[l]->sumV/mlmc_est[l]->N - 
						(mlmc_est[l]->sumY/mlmc_est[l]->N)*(mlmc_est[l]->sumY/mlmc_est[l]->N);
					hl = 1.0/pow(M,l);
					sum += sqrt(mlmc_est[l]->V/hl);
				}
				
				//define optimal Nl
				for (l=0; l<=L; l++) {
					mlmc_est[l]->newN = ceil(2.0 * sqrt(mlmc_est[l]->V * hl) * sum / (eps * eps));
					printf("newN[%d]:%d\n", l, mlmc_est[l]->newN);
				}

				//update sample sum
				for (l=0; l<=L; l++) {
					dNl = mlmc_est[l]->newN - mlmc_est[l]->N;
					printf("dN[%d]:%d\n", l, dNl);
					if (dNl > 0) {
						mlmc_est[l]->N += dNl;
						mlmc_level_l_estimation(pricer, mlmc_est[l], l, M, dNl);
					}
				}

				//test for convergence
				if (extraporation_flag) {
					if (L > 1) {
						con = mlmc_est[L]->sumY/mlmc_est[L]->N - (1.0/(double)n) * mlmc_est[L-1]->sumY/mlmc_est[L-1]->N;
						if (fabs(con) < 1.0*(M*M - 1.0)*eps/sqrt(2.0)) {
							converged_flag = 1;
						}
					}
				} else {
					if (L > 1) {
						con = fmax(fabs(mlmc_est[L]->sumY/mlmc_est[L]->N), fabs(mlmc_est[L-1]->sumY/mlmc_est[L-1]->N/M));
						if (con < (M-1)*eps/sqrt(2.0)) {
							converged_flag = 1;
						}
					}
				}
				printf("N:%d, Y:%lf, V:%12.12lf, sumV:%lf\n", mlmc_est[L]->N, mlmc_est[L]->sumY/mlmc_est[L]->N, mlmc_est[L]->V, mlmc_est[L]->sumV);
				printf("L:%d, con:%lf, test:%lf \n\n", L, con, (M-1)*eps/sqrt(2.0));
			}
			sum=0.0;
			for (l=0; l<=L; l++) {
				sum += mlmc_est[l]->sumY / mlmc_est[l]->N;
			}
			if (extraporation_flag) {
				sum += mlmc_est[L]->sumY/mlmc_est[L]->N/(M-1);
			}
			printf("P:%.12e\n\n", sum);
		}

		for (l=0; l<L_MAX; l++) {
			free(mlmc_est[l]);
		}
		free_SDE_WA_PRICER(pricer);
		free_SDE_WA_SLTN(sl);
		free_SDE_WA_SYSTEM(sde);
	}

	return SDE_WA_SUCCESS;
}


/**
 * MLMC_L allocation
 */
MLMC_L *alloc_level_l_estimator(void) {
	MLMC_L *mlmc_l;

	mlmc_l = (MLMC_L *)malloc(sizeof(MLMC_L));
	mlmc_l->sumY = 0.0;
	mlmc_l->sumV = 0.0;
	mlmc_l->V = 0.0;
	mlmc_l->N = 10000;
	mlmc_l->newN = 0;

	return mlmc_l;
}

void init_level_l_estimator(MLMC_L *mlmc_l) {
	mlmc_l->sumY = 0.0;
	mlmc_l->sumV = 0.0;
	mlmc_l->V = 0.0;
	mlmc_l->N = 10000;
	mlmc_l->newN = 0;
}

void gen_unif_rand_number_gsl_pseudo(void *gen, double *rand, int n) {
	gsl_rng *rng = gen;
	int j=0;

	for (j=0; j<n; j++) {
		rand[j] = gsl_rng_uniform(rng);
	}
}

void gen_unif_rand_number_gsl_quasi(void *gen, double *rand, int n) {
	gsl_qrng *q = gen;

	gsl_qrng_get(q, rand);
}

/**
 * transform uniform random number to normal random number.
 * transform uniform random number to normal random number.
 * @param unif_rand n dimension uniform random number.
 * @param normal_rand n dimension normal random number.
 * @param n dimension of randm number. n must be even number.
 * @return transformed noraml random number
 */
void trans_rand_unif_to_normal(double *unif_rand, double *normal_rand, int n) {
	int i;
	double *u_seq1, *u_seq2;
	double *n_seq1, *n_seq2;

	if ((n % 2) != 0) {
		printf("%d error\n", n);
		return;
	}
	u_seq1 = unif_rand;
	u_seq2 = u_seq1 + n/2;
	n_seq1 = normal_rand;
	n_seq2 = n_seq1 + n/2;
	for (i=0; i<n/2; i++) {
		n_seq1[i] = sqrt(-2.0*log(u_seq1[i]))*cos(2.0*M_PIl*u_seq2[i]);
		n_seq2[i] = sqrt(-2.0*log(u_seq1[i]))*sin(2.0*M_PIl*u_seq2[i]);
	}
}

/**
 * level l estimaton
 * level l esitmator for MLMC
 * @param sl 
 * @param sde
 * @param l level
 * @param M base of number of partitions 
 * @param N num of samles for level l estimation
 * @return estimated expectation
 */
int mlmc_level_l_estimation(SDE_WA_PRICER *pricer, MLMC_L *estimator, int l, int M, int N) {
	SDE_WA_SLTN *sl = pricer->sltn;
	SDE_WA_SYSTEM *sde = sl->sde;
	int i,j,k,m;
	double sum=0.0;
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
		for (sum=0.0, i=0; i < N; i++){
			for (j=0; j<sde->dim_y; j++) {
				xf[0][j] = pricer->init_y[j];
			}
			pricer->next_rand(pricer->gen_unif_rand, sp2, sde->dim_BM*2);
			trans_rand_unif_to_normal(sp2, sp, sde->dim_BM*2);

			next_SDE_WA(sl, pricer->T, xf[0], xf[1], sp);
			tmp_pt=xf[0];
			xf[0]=xf[1];
			xf[1]=tmp_pt;

			estimator->sumY += pricer->payoff_function(xf[0], sde->params);
			estimator->sumV += pricer->payoff_function(xf[0], sde->params) * pricer->payoff_function(xf[0], sde->params);
		} /* for i */
	} else {
		sp=(double *)malloc(sizeof(double)*sde->dim_BM);
		sp2=(double *)malloc(sizeof(double)*sde->dim_BM);
		dt1=pricer->T/(double)n;
		dt2=M/(double)n;
		for (sum=0.0, i=0; i < N; i++){
			pricer->next_rand(pricer->gen_unif_rand, u_seq, n*sde->dim_BM);
			trans_rand_unif_to_normal(u_seq, n_seq, n*sde->dim_BM);

			for (j=0; j<sde->dim_y; j++) {
				xf[0][j] = pricer->init_y[j];
				xc[0][j] = pricer->init_y[j];
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
			tmp = pricer->payoff_function(xf[0], sde->params) - pricer->payoff_function(xc[0], sde->params);
			estimator->sumY += tmp; 
			estimator->sumV += tmp*tmp;
		} /* for i */
	}
	return sum;
}

