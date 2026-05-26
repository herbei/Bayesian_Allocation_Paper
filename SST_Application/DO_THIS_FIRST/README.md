# SST Data Download Quick Start

Run these steps before rebuilding the SST application inputs, the ocean graph, MCMC runs, or figures.

All scripts in this folder resolve paths against the parent `Rcode/SST_Application` directory, so they can be run either from that parent directory or from `DO_THIS_FIRST`.

## Required R Packages

```r
install.packages(c("ncdf4", "Matrix", "maps"))
```

## 1. Download Observations and Base Grid Inputs

From `Rcode/SST_Application`:

```sh
Rscript DO_THIS_FIRST/data_retrieval_sst.R
```

This downloads:

- NOAA OISST monthly observations to `Data/raw/oisst_v2_sst_mnmean.nc`

It writes the observation and grid files to `Data/processed/`, including `pacific_Yo.csv`, `pacific_sst_grid.csv`, `pacific_land_mask.csv`, and `pacific_sst_data.rds`.

## 2. Download E3SM Ensemble Members

For the 21-member set used by the current processed inputs:

```sh
Rscript DO_THIS_FIRST/download_e3sm_ensembles.R
```

Raw ensemble files are written to `Data/raw/tos_Omon_*.nc`. The default member set is `r1i1p1f1` through `r21i1p1f1`. Use `--help` to see options for alternate members, versions, output directories, or time ranges.

## 3. Extract Ensemble Inputs

After the raw ensemble NetCDF files are present:

```sh
Rscript DO_THIS_FIRST/data_retrieval_ens.R
```

This writes the ensemble inputs under `Data/processed/`, especially `pacific_Yc_ensemble_members/pacific_Yc_<member>.csv`, which is what the ensemble MCMC workflow reads. It also writes `pacific_Yc.csv` as the ensemble mean for scripts that expect a single computer-output vector.

## 4. Build the Ocean Graph and Laplacian

After `pacific_sst_grid.csv` and `pacific_land_mask.csv` exist:

```sh
Rscript DO_THIS_FIRST/compute_ocean_laplacian.R --force-full-eigen
```

This writes `pacific_ocean_graph_laplacian.rds`, `pacific_ocean_edges.csv`, `pacific_ocean_sites_with_local_index.csv`, and `pacific_ocean_laplacian_eigendecomp.rds` under `Data/processed/`.

The full eigendecomposition is large and needed by the MCMC workflow. For a graph-only rebuild, use:

```sh
Rscript DO_THIS_FIRST/compute_ocean_laplacian.R --no-full-eigen
```

## Source Details

See `DATA_SOURCES_SST.md` in this folder for data hosts, URLs, variables, grid conventions, and processing details.

