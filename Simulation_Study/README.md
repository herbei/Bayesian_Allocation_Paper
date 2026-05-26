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

Generated figure outputs are intentionally excluded.

## Dependencies

The MCMC entry point uses base R plus the fixed files in `Data/inputs/`.

The figure scripts require these R packages:

```r
install.packages(c("ggplot2", "gridExtra", "maps", "scales"))
```

PDF outputs are cropped with the external `pdfcrop` command, which is usually
provided by TeX Live/MacTeX. Check it with:

```sh
command -v pdfcrop
```

## Run

From this directory:

```sh
Rscript run_mcmc.R
```

By default this writes `Data/results/simulation_study_selected_run.RData`.

To make a short smoke-test run:

```sh
N_ITER=10 BURN_IN=2 MCSTORE=1 Rscript run_mcmc.R
```

Using the selected MCMC artifact, recreate the Section 4 figures:

```sh
./run_section4_figures.sh Data/results/simulation_study_selected_run.RData
```

The main environment-variable overrides are `N_ITER`, `BURN_IN`, `MCSTORE`,
`PROGRESS_EVERY`, `SEED`, `K0`, `M_ENS`, `P_K`, `DELTA1`, `DELTA2`,
`DELTA01`, `DELTA02`, `SIG2C_INIT`, `SIG2O_INIT`, `UPDATE_SIGMA2_AFTER`,
and `OUTPUT_FILE`.
