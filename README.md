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

The installed Calibr revision `c01a645d1f2b1cd142c079008a71ad9f00d82c56`
has a GLMM group-name bug: it incorrectly labels group-specific estimates as
global. For this revision, apply the included patch before running GLMMs.
The following PowerShell commands install a corrected copy into a separate R
library without changing a shared installation:

```powershell
$src = Join-Path $env:TEMP ("Calibr-" + [guid]::NewGuid())
$env:R_LIBS_USER = Join-Path $env:USERPROFILE "R-calibr-fixed"
New-Item -ItemType Directory -Force $env:R_LIBS_USER | Out-Null
git clone https://github.com/marcnadon/Calibr.git $src
git -C $src checkout c01a645d1f2b1cd142c079008a71ad9f00d82c56
git -C $src apply (Resolve-Path "patches/calibr-glmm-group-names.patch").Path
R.exe CMD INSTALL --library=$env:R_LIBS_USER $src
Rscript tests/test_calibr_glmm.R
```

GLMM also requires a TMB version compatible with the installed `glmmTMB`;
resolve any package version-mismatch warning before running. The patched
Calibr library must be on `.libPaths()` for every calibration invocation.
Existing GLMM exports created with the affected revision should be rerun.

Run scripts from the repository root in this order:

```powershell
Rscript script/0_compare_nspc_data.R
Rscript script/1_run_calibr.R
Rscript script/2_coefficients.R
Rscript script/3_plots.R
```

`1_run_calibr.R` uses I-only nSPC data for calibration and runs both GLM and
GLMM for abundance and biomass across the five regions and belt/tow methods.
Pass a species code to refresh only that group while retaining
other groups in the calibration exports. Optional metric, region, model, and
method arguments can narrow the run further.
Set the `CALIBR_OUTPUT_DIR` environment variable to an existing directory to
write calibration exports somewhere other than `output`. Species-only runs
require existing exports in the selected directory.
The current all-group run also provides `output/calibration_status.csv`:
`group_specific` means a finite group-specific coefficient was exported,
`global_fallback` means the GLMM used its global coefficient,
`no_secondary_summary` means no secondary-method result was exported, and
`invalid_coefficient` identifies a missing, non-finite or non-positive
positive-density coefficient. Calibr rounds GLMM coefficients to two decimal
places, so a small positive coefficient can appear as zero in the export and
must not be used as a divisor.
`2_coefficients.R` summarizes the GLM results and accepts optional species
codes to limit which coefficient plots are refreshed. `3_plots.R` uses the
I-only reference so its calibrated outputs match the calibration coefficients;
it defaults to abundance/APVI/SAMOA and accepts an optional metric and region:

```powershell
Rscript script/1_run_calibr.R APVI
Rscript script/2_coefficients.R APVI
Rscript script/1_run_calibr.R APVI abund MARIAN GLMM tow
Rscript script/3_plots.R abund NWHI
Rscript script/3_plots.R biom NWHI
Rscript script/3_plots.R abund NWHI maps
Rscript script/3_plots.R biom MARIAN maps CAIG
```

The optional `maps` selection refreshes only the four map outputs using the
existing calibrated dataset; a fourth argument selects a species instead of
the default APVI. Run the full plotting command first if that dataset is
missing or needs updating.

Maps use the R `sf` package and the bundled Natural Earth 1:10m land polygons
([Natural Earth](https://www.naturalearthdata.com/), public domain) for
regional views. Island-detail panels use 2025 EOxCloudless imagery (modified
Copernicus Sentinel data) from [EOX](https://cloudless.eox.at/) under
[CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/).
The imagery and maps containing it are for non-commercial use and subject to
share-alike requirements; check the [EOX terms](https://cloudless.eox.at/license-non-commercial)
before reuse. Imagery is historical, may have gaps, and is not navigation-grade.
These panels require the R `curl`, `jpeg` and `png` packages and an internet connection when the
imagery is not already cached. Whole-degree public map tiles are downloaded
and cached locally in `data/.eox_2025_cache/`; precise survey coordinates are
not sent to the imagery service. If EOX reports that a high-resolution tile is
unavailable, a lower-resolution tile is used with a warning; other download
errors or invalid tiles stop the map run.
The composite map (`map_b`) combines the regional overview with island-detail
panels, each with its own scale. The overview shows means of positive
observations in 0.1-degree cells; its legend is separate from the island legend.
The island-detail map averages all composite observations at each location,
including zeros, and shows zero means as white points. Additional close-ups
supplement full-island views when outlying coordinates compress the map.
The year-by-method map (`map_d`) groups years by survey method, with a shared
regional extent and legend. It shows mean positive observations in 0.1-degree
cells; use the island-detail panels to inspect small or widely separated islands.

The comparison script uses IUN for the standard time series.

`0_prep_nspc_data.R` is the older raw-data preparation workflow. Its generated
values are marked provisional because historical normalization remains
unresolved; do not use those outputs as a substitute for the supplied
site-summary files.
