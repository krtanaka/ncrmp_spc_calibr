# SPC calibration

NCRMP survey-data calibration using the R Calibr package.

## Updated nSPC summaries

Place the supplied Parquet files in `data/tanaka_2026/`. The workflow uses:

- `site_summary_SPECIES_I.parquet` as the nSPC calibration reference. This
  contains instantaneous observations only.
- `site_summary_SPECIES_IUN.parquet` for standard NCRMP time-series summaries.
  IUNFT is intentionally excluded from calibration because F/T observations
  represent additional effort and are not directly comparable.

The summaries already apply the site-level species grouping, species-code
combination, cryptic/non-fish exclusions, replicate averaging, and biomass
length-conversion factor. The workflow uses `DENSITY_M2` and `BIOMASS_G_M2`;
it does not rebuild these site values from observation-level records.
For 2026 rows, the supplied `DATE_` timestamp carries year `0026`; the reader
repairs that year from the matching `OBS_YEAR` value before calibration and
plotting.

Install the R `arrow` package to read the Parquet summaries:

```r
install.packages("arrow")
```

If Calibr is not installed, install the upstream package with
`remotes::install_github("marcnadon/Calibr")`.

Run scripts from the repository root in this order:

```powershell
Rscript script/0_compare_nspc_data.R
Rscript script/1_run_calibr.R
Rscript script/2_coefficients.R
Rscript script/3_plots.R
```

`1_run_calibr.R` uses I-only nSPC data for calibration and runs both GLM and
GLMM for abundance and biomass across the five regions and belt/tow methods.
`2_coefficients.R` summarizes the GLM results. `3_plots.R` uses the I-only
reference so its calibrated outputs match the calibration coefficients; its
current selection is abundance/APVI/SAMOA. The comparison script uses IUN for
the standard time series.

`0_prep_nspc_data.R` is the older raw-data preparation workflow. Its generated
values are marked provisional because historical normalization remains
unresolved; do not use those outputs as a substitute for the supplied
site-summary files.
