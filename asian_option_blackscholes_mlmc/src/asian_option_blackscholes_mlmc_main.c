/**
 * @file    asian_option_blackscholes_mlmc_main.c
 * @brief   Main function of MLMC Method to calculate asian options under blackscholes model.

 *          Main function of MLMC Method to calculate asian options under blackscholes model. 
 * @date    09/Mar/2015 (Mon) 17:32 
 * @author  i05nagai
 */
#include <stdlib.h> /* for malloc */
#include <stdio.h>
#include <math.h>
#include <sde_wa/sde_wa.h>
#include <sde_wa/sde_wa_em.h>
#include "asian_option_blackscholes.h"
#include "asian_option_blackscholes_mlmc.h"
#include <sde_mlmc/mlmc.h>
#include <mt19937/mt19937.h>


int main(void){
	//enum ALG alg = E_M;
	//enum ALG alg = N_V;
	enum ALG alg = N_N;
	int m_moment = 7; /* Romberg */
	int M = 4;
	struct ABS_params abs_params;
	double mu=0.05;
	double sigma=0.2;
	double T=1.0;
	double K=1.0;
	double x0=1.0;
	double x1=0.0;
	abs_params.mu=mu;
	abs_params.sigma=sigma;
	abs_params.T=T;
	abs_params.K=K;
	MLMC_L *mlmc_ests;
	{
		unsigned long init[4]={0x123, 0x234, 0x345, 0x456};
		unsigned long length=4;
		init_by_array(init, length);

		SDE_WA_SYSTEM *sde;
		SDE_WA_SLTN *sl;
		SDE_MLMC *mlmc;
		double eps=0.0001;
		int l;
		int L;
		sde=alloc_SDE_WA_SYSTEM(2, 1, &abs_params);
		sde->sde_type=STR;

		sde->V[0]=abs_str_V_0;
		sde->V[1]=abs_V_1;

		sde->drift_corrector[0]=NULL;
		sde->drift_corrector[1]=diff_abs_V_1; /* for EM */

		sde->exp_sV[0]=exp_abs_sVy_0; /* for NV */
		sde->exp_sV[1]=exp_abs_sVy_1; /* for NV */

		sl=alloc_SDE_WA_SLTN(alg, m_moment, sde);

		mlmc=alloc_SDE_MLMC(sl);
		mlmc->payoff_function=asian_option_blackscholes_call_payoff;
		if (sl->alg == E_M) {
			printf("Discretization Mehtod:E_M\n");
			mlmc->level_l_estimation=mlmc_level_l_estimation_em;
		} else if(sl->alg == N_N) {
			printf("Discretization Mehtod:N_N\n");
			mlmc->level_l_estimation=mlmc_level_l_estimation_nn;
		} else {
			printf("Discretization Mehtod:N_V\n");
			mlmc->level_l_estimation=mlmc_level_l_estimation_nv;
		}
		mlmc->init_y[0] = x0;
		mlmc->init_y[1] = x1;
		mlmc->T = 1.0;

		mlmc_ests = mlmc->list;
		for (l=0; l<=-4; l++) {
			init_level_l_estimator(mlmc_ests);
			mlmc_ests->N = 40000;
			mlmc->level_l_estimation(mlmc, mlmc_ests, l, M, mlmc_ests->N);
			mlmc_ests->V = mlmc_ests->sumYY/mlmc_ests->N - (mlmc_ests->sumY/mlmc_ests->N)*(mlmc_ests->sumY/mlmc_ests->N);
			printf("L:%d M:%d\n", l, M);
			printf("N:%d, Y:%12e, V:%12.12lf, sumYY:%lf\n", mlmc_ests->N, mlmc_ests->sumY/mlmc_ests->N, mlmc_ests->V, mlmc_ests->sumYY);
			if (mlmc_ests->next == NULL) {
				add_level_l_estimator(mlmc_ests);
			}
			mlmc_ests = mlmc_ests->next;
		}
		//exit(-1);

		int i;
		double epss[] = {0.001, 0.0005, 0.0002, 0.0001, 0.00005};
		for (i=0; i<5; i++) {
			L=-1;
			eps = epss[i];
			mlmc_ests = mlmc->list;
			while (mlmc_ests != NULL) {
				init_level_l_estimator(mlmc_ests);
				
				mlmc_ests = mlmc_ests->next;
			}
			printf("eps:%lf\n", eps);

			//MLMC Algorithm
			mlmc_algorithm(mlmc, M, eps, 0);

			mlmc_ests = mlmc->list;
			while (mlmc_ests->next != NULL) {
				printf("N:%d, V:%12e\n", mlmc_ests->N, mlmc_ests->V);
				mlmc_ests = mlmc_ests->next;
			}
			printf("P:%lf\n", mlmc->P);
		}

		free_SDE_MLMC(mlmc);
		free_SDE_WA_SLTN(sl);
		free_SDE_WA_SYSTEM(sde);
	}

	return SDE_WA_SUCCESS;
}


