SDE\_WA\_MLMC Samples
====
## Overview

This repository consists of some samples of the library([SDE\_WA\_MLMC](https://github.com/i05nagai/SDEWeakApproximation)).

## Description

The samples estimate option pricing problems under two models discretizing by Kusuoka schemes and Euler Maruyama scheme.
Options, models and discretization schemes  are as follows.

Two options:

* European Options
* Asian Options

Two models:

* Black-Scholes Model
* Heston Model

Discretization schemes:

* Kusuoka Schemes
    * [Ninomiya-Victoir](http://www.tandfonline.com/doi/abs/10.1080/13504860701413958#.VNm8ZsbAJDQ)
	* [Ninomiya-Ninomiya](http://link.springer.com/article/10.1007/s00780-009-0101-4)
* Stochastic Taylor Expansions
    * Euler Maruyama Scheme

## Requirement

This samples depend on a library([SDE\_WA\_MLMC](https://github.com/i05nagai/SDEWeakApproximation)).

In MacOSX and Linux, you need to make a path or create a symbolic link to `libsde_wa.a`, which is the library of SDE\_WA\_MLMC, as below.

	`ln -s SDEWeakApproximation/libsde_wa.a SDEWeakApproximationMLMCSamples/libsde_wa.a`

## Usage

1. Move to a sample directory(i.e. asian\_option\_heston\_mlmc)

    `cd asian_option_heston_mlmc`

2. Create Makefile

    `make -f basic.mk Makefile`

3. Make an execution(i.e. asian\_option\_heston\_mlmc\_main)

    `make`

4. Execute

    `./asian_option_heston_mlmc_main`

## Licence

This software is released under the MIT License, see LICENSE.

## Author

[i05nagai](https://github.com/i05nagai)

