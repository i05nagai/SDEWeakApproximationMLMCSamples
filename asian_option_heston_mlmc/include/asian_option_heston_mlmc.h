/**
 * @file    asian_option_heston_mlmc.h
 * @brief   header file of functions of level l estimation.

 *          header file of functions of level l estimation. 
 * @date    09/Mar/2015 (Mon) 17:32 
 * @author  i05nagai
 */
#ifndef _SDE_AH_MLMC_
#define _SDE_AH_MLMC_
#include <sde_mlmc/mlmc.h>

int mlmc_level_l_estimation_em(SDE_MLMC *mlmc, MLMC_L *estimator, int l, int M, int N);
int mlmc_level_l_estimation_nv(SDE_MLMC *mlmc, MLMC_L *estimator, int l, int M, int N);
int mlmc_level_l_estimation_nn(SDE_MLMC *mlmc, MLMC_L *estimator, int l, int M, int N);

#endif

