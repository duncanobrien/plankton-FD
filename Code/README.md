# plankton-FD Code

## Recommended code workflow

[*calc_system_states.R*](calc_system_states.R) - from the raw plankton abundances/densities, estimate ecosystem state using five separate measures: total abundance, community principal component, Fisher information, multivariate variance index and trophic ratio.

[*calc_FD.R*](calc_FD.R) - from the raw plankton abundances/densities and trait information, estimate functional diversity (FD).

[*calc_permuted_ccf.R*](calc_permuted_ccf.R) - estimate permuted cross correlations between the detrended FD and system state measures.

[*calc_permuted_ccm.R*](calc_permuted_ccm.R) - estimate permuted convergent cross mapp skills between the detrended FD and system state measures.

[*Figure1.R*](Figure1.R), [*Figure2.R*](Figure2.R), [*Figure3.R*](Figure3.R), [*Figure4.R*](Figure4.R), [*supplementary_figures.R*](supplementary_figures.R) - generate figures.

## Supporting functions
[*fisher_information*](fisher_information) - collection of functions to estimate Fisher Information. See subfolder for details.

[*mvi_fn.R*](mvi_fn.R) - estimates the multivariate variance index following Brock, W. A., and S. R. Carpenter. 2006. Variance as a leading indicator of regime shift in ecosystem services. *Ecology and Society* 11(2): 9. http://www.ecologyandsociety.org/vol11/iss2/art9/.

[*tidy_FD_fn.R*](tidy_FD_fn.R) - estimates functional diversity using the excellent mFD package (https://doi.org/10.1111/ecog.05904) and [*melodic_fn.R*](melodic_fn.R).

[*melodic_fn.R*](melodic_fn.R) - melodic function published in de Bello, F., Carmona, C.P., Lepš, J. et al. 2016. Functional diversity through the mean trait dissimilarity: resolving shortcomings with existing paradigms and algorithms. *Oecologia* 180, 933–940. https://doi.org/10.1007/s00442-016-3546-0 for abundance-weighted MPD and Rao's Q.

[*diff_perm_ccf_fn.R*](diff_perm_ccf_fn.R) - permutation cross correlation function allowing various detrending methods, permutation methods, and lag ranges to be specified.

[*ccm_perm_fn.R*](ccm_perm_fn.R) - permutation convergent cross mapping function using the rEDM package (https://CRAN.R-project.org/package=rEDM) allowing various detrending methods and lag ranges to be specified.
