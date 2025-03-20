
# Stokes-NNKIFMM
[![License: AGPL v3](https://img.shields.io/badge/License-AGPL_v3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)

Stokes-NNKIFMM is an adaptive particle kernel-independent FMM Code based on [^1] for the fast computation of Regulasied Stokeslet dynamics [^2]. 
The Code employs the nearest neighbour mapping (detailed in [^3][^4]) and Richardson extrapolation (detailed in [^5]) to quickly evaluate the Mobility equations for microscale biological fluid dynamics.
Mobility equations are solved using GMRES [^6] using a preconditioner inspired by [^7].
This work was done as part of my master of physics at the University of Birmingham supervised by David Smith and supported by Meurig Gallagher. This work has not been published but my Dissertation can be found [Here](https://www.github.com/Walkersb101/Dissertation)

## How to Use
The Code supports Serial, multi-threading and GPU (including multiple GPUs) on a single node. 
Test cases can be found in the [Tests](Tests) folder, The scripts used to generate data used in the figures in my Dissertation can be found in the [Data Generation](Data%20Generation) folder and the code used to plot this data can be found in the [Figure Plotting](Figure%20Plotting) folder.

## References
[^1]: Ying, L., Biros, G., & Zorin, D. (2004). A kernel-independent adaptive fast multipole algorithm in two and three dimensions. Journal of Computational Physics, 196(2), 591-626.
[^2]: R. Cortez and S. J. Sci Comput. “The Method of Regularized Stokeslets”. In: Artic. SIAM J. Sci. Comput. 23.4 (2001), pp. 1204–1225.
[^3]: Gallagher MT, Smith DJ. Meshfree and efficient modelling of swimming cells. Physical Review Fluids. 2018 May 31;3(5):053101.
[^4]: Smith DJ. A nearest-neighbour discretisation of the regularized stokeslet boundary integral equation. Journal of Computational Physics. 2018 Apr 1;358:88-102.
[^5]: Gallagher M. T. and Smith D. J. 2021 The art of coarse Stokes: Richardson extrapolation improves the accuracy and efficiency of the method of regularized stokesletsR. Soc. Open Sci.8210108
[^6]: Y. Saad and M. H. Schultz. “GMRES: A Generalized Minimal Residual Algorithm for Solving Nonsymmetric Linear Systems”. In: SIAM J. Sci. Stat.Comput. 7.3 (1986), pp. 856–869.
[^7]: M. W. Rostami and S. D. Olson. “Fast algorithms for large dense matrices with applications to biofluids”. In: J. Comput. Phys. 394 (Oct. 2019), pp. 364–384.
