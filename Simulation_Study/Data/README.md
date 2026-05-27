# Simulation Study Data

`inputs/` contains the fixed data files needed to reproduce the Section 4
simulation study:

- `phi_o_data.txt`: observed OXY field on the full grid before subsetting to
  observation locations.
- `api.txt`: one-based grid indices for the observation locations.
- `oxy_ensemble_20_members.rds`: 20-member computer-model ensemble.
- `oxy_full_grid_laplacian_eigendecomp.rds`: cached eigendecomposition of the
  19 x 37 full-grid graph Laplacian.
- `trueXX.txt`: true latent field used by the figure scripts and MCMC-run
  comparison summaries.
- `Binit.txt`: selected initial allocation vector for the reported run.

`results/` contains the selected MCMC artifact,
`simulation_study_selected_run.RData`, and is also where new MCMC artifacts and
generated downstream summaries are written by default. Running
`../compare_mcmc_runs.R` writes comparison outputs to
`results/comparison/` for the folder selected by `mcmc_runs_dir`.
