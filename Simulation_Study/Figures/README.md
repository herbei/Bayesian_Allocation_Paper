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

Run the scripts from `src/`; paths are resolved relative to the submitted code
directory. By default, the scripts write generated figures to `out/`. Set
`FK_OUTPUT_DIR` to write Figures 01--03 elsewhere, and set
`OXY_FIGURE_CHAIN_FILE` or the figure-specific chain-file variables for
Figures 04--06.

```sh
Rscript src/Figure01.R
Rscript src/Figure02.R
Rscript src/Figure03.R
Rscript src/Figure04.R
Rscript src/Figure05.R
Rscript src/Figure06.R
```

All fixed inputs are read from `../Data/inputs/`.
