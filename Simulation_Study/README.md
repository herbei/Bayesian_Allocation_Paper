# Section 4: Simulation Study

This folder contains the code used for the oxygen simulation study in Section 4.

## Contents

- `run_mcmc.R`: selected-run MCMC entry point matching the reported compromise run
  (`K = 703`, `K0 = 300`, `p_k = 0.5`, `delta1 = delta2 = 2.5`).
- `modules/spectral2_ens.R`: spectral modal-allocation MCMC implementation.
- `run_section4_figures.sh`: recreates the Section 4 figure panels from an MCMC artifact.
- `Data/inputs/`: fixed inputs used by the MCMC, simulator, and figure scripts.
- `Data/results/`: MCMC artifacts and downstream generated summaries.
- `Figures/src/`: manuscript figure scripts and simulator helpers.
- `Figures/out/`: generated figure outputs.

Generated outputs are intentionally excluded.

## Run

From this directory:

```sh
Rscript run_mcmc.R
```

To make a short smoke-test run:

```sh
N_ITER=10 BURN_IN=2 MCSTORE=1 Rscript run_mcmc.R
```

After an MCMC artifact exists, recreate the Section 4 figures:

```sh
./run_section4_figures.sh Data/results/simulation_study_selected_run.RData
```

The main environment-variable overrides are `N_ITER`, `BURN_IN`, `MCSTORE`,
`PROGRESS_EVERY`, `SEED`, `K0`, `P_K`, `DELTA1`, `DELTA2`, `DELTA01`, and `DELTA02`.
