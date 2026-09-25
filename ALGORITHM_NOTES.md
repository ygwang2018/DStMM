# DStMM algorithm-to-code map

The implementation in `R/dstmm.R` follows the manuscript equations directly.

| Manuscript step | Code |
|---|---|
| Path recursion for `mu_s`, `alpha_s`, `Sigma_s` (13–15) | `.dstmm_collapse_path()` |
| Exact GHST pathway density (19–21) | `.dstmm_logdghst()`, `.dstmm_estep()` |
| Path responsibilities `tau_js` (23) | `.dstmm_estep()` |
| GIG law of `H | y,s` (28) | `.dstmm_rgig()`, `.dstmm_draw_latent_one()` |
| Student-t/IG limit (29–30) | `.dstmm_logdmvt()`, IG branch in E-step/draw |
| GIG moments (34–35) | `.dstmm_gig_moments()` |
| Conditional latent-state simulation (36–38) | `.dstmm_draw_latent_one()` |
| Weighted regression update (39–44) | `.dstmm_mstep()` |
| Layer mixing proportions (45) | `.dstmm_mstep()` |
| Diagonal `Psi` update (46) | `.dstmm_mstep()` |
| Degrees-of-freedom root equation (47) | `.dstmm_update_nu()` |
| Complete stochastic/MCEM loop (Algorithm 1) | `dstmm()` |

## Parameter naming

The uploaded RDMM code uses `H` for factor-loading matrices and a Gamma precision for its robust scale representation. To avoid collision with the DStMM manuscript's inverse-gamma variable `H`, the new implementation uses:

- `lambda` for loading matrices `Lambda`;
- `eta` for local intercepts;
- `delta` for local skewness shifts;
- `psi` for specific covariance matrices;
- `H` only for the observation-path inverse-gamma scale inside the DStMM sampler.

## Nested-model checks

With every `delta` set to zero, the DStMM pathway density switches to the multivariate Student-t limit. Thus the code preserves the theoretical DStMM → RDMM nesting at the observed-path density level. The fitting scripts explicitly include `kappa = 0` so the extra skewness parameters can be checked through their estimated RMSE and clustering recovery.
