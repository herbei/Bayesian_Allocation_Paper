# Simulation Study Figures

This folder contains the figure scripts and shared plotting helpers used to
reproduce the Section 4 figures shipped in
`Spatial_Statistics/manuscript/figures/`.

- `src/Figure01.R`: reproduces manuscript `Figure01.pdf`, four representative
  computer-model ensemble members.
- `src/Figure02.R`: reproduces manuscript `Figure02.pdf`, the pointwise ensemble
  variance map and variance histogram.
- `src/Figure03.R`: reproduces manuscript `Figure03.pdf`, the true latent field and
  induced discrepancy field.
- `src/Figure04.R`: reproduces manuscript `Figure04.pdf`, posterior means for
  the latent field and model bias.
- `src/Figure05.R`: reproduces manuscript `Figure05.pdf`, posterior standard
  deviations.
- `src/Figure06.R`: reproduces manuscript `Figure06.pdf`, posterior allocation
  probabilities.
- `src/plot_utils.R`: shared OXY data-loading, coastline, plotting, PDF
  cropping, panel-layout, chain-file, and output-directory helpers used by the
  figure scripts.
- `out/`: default output directory for generated figures.

## Dependencies

The figure scripts require the R packages `ggplot2`, `gridExtra`, `maps`, and
`scales`. PDF cropping requires the external `pdfcrop` command.

## Run

From this `Figures/` directory:

```sh
Rscript src/Figure01.R
Rscript src/Figure02.R
Rscript src/Figure03.R
Rscript src/Figure04.R
Rscript src/Figure05.R
Rscript src/Figure06.R
```

The scripts resolve paths relative to the submitted code directory, so they can
also be run through `../run_section4_figures.sh`. By default, Figures 01--03
write to `out/`; Figures 04--06 write to `out/<chain-file-name>/` because they
depend on an MCMC artifact.

Set `FK_OUTPUT_DIR` to write Figures 01--03 elsewhere. Set
`OXY_FIGURE_CHAIN_FILE` or the figure-specific chain-file variables
`OXY_FIGURE04_CHAIN_FILE`, `OXY_FIGURE05_CHAIN_FILE`, and
`OXY_FIGURE06_CHAIN_FILE` for Figures 04--06. Set `OXY_FIGURE_OUTPUT_DIR` to
override the posterior-figure output directory.

All fixed inputs are read from `../Data/inputs/` by default.
