# dofcalcs 
This package provides some functions for calculating degrees of freedom for
models with with glmmTMB. It is constructed for use with downstream package 
[`powerutilities`](https://www.github.com/SimonShamusRiley/powerutilities).
All the code is adapted from work by Josh Price which
was first published as part of a `glmmTMBbasicDFs` repository and later
incorporated into the [lmmaov](https://github.com/jprice80/lmmaov) repository.

In keeping with the requirements of the AGPL-3 licence under which the code
was originally published, this package is likewise published under the AGPL-3
licence, and every effort will be made here to provide details of the
modifications made to the original code:

- The functions `residual_aov` has been added, which assigns DFs based on the
  value returned by `df.residual`.

- The function `asymptotic_aov` has been added, which retains the use of 
  infinite degrees of freedom as implemented in `car::Anova`

- The `satterthwaite_aov` function has been added, which calculates using
  DFs using the `glmmTMB::dof_satt` function.

- The `anova`, `inner_outer_aov` and `inner_outer_summary` functions were 
  removed (an `anova` method for `glmmTMB` objects is instead implemented in
  the powerutilities package).