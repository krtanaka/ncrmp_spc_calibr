# ==============================================================================
# script/0_prep_nspc_data.R
#
# Process raw survey data from fish-paste (through 2024) and generate regional
# nSPC site-level abundance and biomass RDS files preserving exact schema,
# column types, and factor levels for downstream calibration models.
# Run from the repository root. Historical normalization remains unresolved:
# generated densities are provisional and must not yet be used for calibration.
# ==============================================================================

library(data.table)

check_rows <- function(dat, invalid, label) {
  if (any(invalid)) {
    stop(label, ": ", sum(invalid), " invalid rows; example SITEVISITID(s): ",
         paste(head(unique(dat$SITEVISITID[invalid]), 5), collapse = ", "),
         call. = FALSE)
  }
}

preserve_factor <- function(values, reference) {
  values <- as.character(values)
  new_levels <- sort(setdiff(unique(values[!is.na(values)]), levels(reference)))
  factor(values, levels = c(levels(reference), new_levels),
         ordered = is.ordered(reference))
}

cat("=== Starting nSPC Data Prep (2009-2024) ===\n")
warning("Historical normalization is unresolved. Outputs are provisional; ",
        "do not use them for calibration yet.", call. = FALSE)

# Path configuration
raw_data_path <- "C:/Users/Kisei.Tanaka/fish-paste/data/ALL_REA_FISH_RAW.rdata"
survey_master_path <- "C:/Users/Kisei.Tanaka/fish-paste/data/SURVEY MASTER.csv"

if (!file.exists(raw_data_path)) {
  stop("Raw fish data file not found at: ", raw_data_path)
}
if (!file.exists(survey_master_path)) {
  stop("SURVEY MASTER file not found at: ", survey_master_path)
}

cat("Loading SURVEY MASTER...\n")
sm <- fread(survey_master_path, select = c("SITEVISITID", "SEC_NAME", "ANALYSIS_YEAR"))
if (anyNA(sm$SITEVISITID) || anyDuplicated(sm$SITEVISITID)) {
  stop("SURVEY MASTER must contain one non-missing row per SITEVISITID.")
}

cat("Loading ALL_REA_FISH_RAW...\n")
load(raw_data_path) # loads `df`
setDT(df)
required <- c("SITEVISITID", "REPLICATEID", "METHOD", "REGION", "OBS_YEAR",
              "EXCLUDE_FLAG", "TRAINING_YN", "OBS_TYPE", "SPECIES", "SIZE_",
              "COUNT", "LW_A", "LW_B", "LENGTH_CONVERSION_FACTOR", "TAXONNAME",
              "ISLAND", "REEF_ZONE", "DEPTH", "DEPTH_BIN", "LATITUDE",
              "LONGITUDE", "DATE_")
missing_cols <- setdiff(required, names(df))
if (length(missing_cols)) {
  stop("Raw data is missing columns: ", paste(missing_cols, collapse = ", "))
}
numeric_cols <- c("SIZE_", "COUNT", "LW_A", "LW_B", "LENGTH_CONVERSION_FACTOR")
if (!all(vapply(df[, ..numeric_cols], is.numeric, logical(1)))) {
  stop("Size, count, and biomass coefficient columns must be numeric.")
}

cat("Filtering valid nSPC survey records (2009-2024)...\n")
# Filter nSPC surveys
nspc <- df[
  METHOD %in% c("nSPC", "nSPC-CCR") &
  OBS_YEAR >= 2009 & OBS_YEAR <= 2024 &
  !is.na(EXCLUDE_FLAG) & EXCLUDE_FLAG == 0 &
  (is.na(TRAINING_YN) | TRAINING_YN == FALSE | TRAINING_YN == 0)
]
rm(df); gc()

# Merge SEC_NAME and ANALYSIS_YEAR from SURVEY MASTER if missing or NA
nspc <- merge(nspc, sm, by = "SITEVISITID", all.x = TRUE, suffixes = c("", ".sm"))
if ("SEC_NAME.sm" %in% names(nspc)) {
  nspc[is.na(SEC_NAME) | SEC_NAME == "", SEC_NAME := SEC_NAME.sm]
  nspc[, SEC_NAME.sm := NULL]
}
if ("ANALYSIS_YEAR.sm" %in% names(nspc)) {
  nspc[is.na(ANALYSIS_YEAR) | ANALYSIS_YEAR == "", ANALYSIS_YEAR := ANALYSIS_YEAR.sm]
  nspc[, ANALYSIS_YEAR.sm := NULL]
}

# Convert to character before assigning a potentially new species code.
nspc[, SPECIES := as.character(SPECIES)]
# Standardize species names: convert CHPA to CHAG
nspc[SPECIES == "CHPA", SPECIES := "CHAG"]

# Calculate 10-cm size bins: (0,10], (10,20], ..., (290,300]
nspc[, SIZE_10cm := as.character(cut(SIZE_, breaks = seq(0, 300, 10)))]

# Filter valid fish observation types for abundance/biomass
# (U: unconstrained, I: instantaneous, N: non-instantaneous)
nspc_obs <- nspc[OBS_TYPE %in% c("U", "I", "N")]

# Target regions
regions <- c("MARIAN", "MHI", "NWHI", "PRIAs", "SAMOA")

# Validate all regions before writing any output.
survey_records <- nspc[REGION %in% regions]
if (length(setdiff(regions, unique(survey_records$REGION)))) {
  stop("Raw data does not contain valid surveys for every target region.")
}
check_rows(survey_records,
           is.na(survey_records$SITEVISITID) | is.na(survey_records$REPLICATEID),
           "Missing survey or cylinder identifiers")
metadata_cols <- c("OBS_YEAR", "ANALYSIS_YEAR", "ISLAND", "SEC_NAME",
                   "REEF_ZONE", "DEPTH", "DEPTH_BIN", "LATITUDE",
                   "LONGITUDE", "DATE_")
for (col in metadata_cols) {
  values <- survey_records[[col]]
  invalid <- is.na(values) | trimws(as.character(values)) == ""
  if (is.numeric(values)) invalid <- invalid | !is.finite(values)
  check_rows(survey_records, invalid, paste("Missing or invalid", col))
}
observations <- nspc_obs[REGION %in% regions]
check_rows(observations,
           is.na(observations$SPECIES) | trimws(observations$SPECIES) == "",
           "Missing species codes")
check_rows(observations,
           is.na(observations$TAXONNAME) |
             trimws(as.character(observations$TAXONNAME)) == "",
           "Missing taxon names")
check_rows(observations, is.na(observations$SIZE_10cm),
           "Missing size or size outside the supported (0,300] cm range")
check_rows(observations,
           !is.finite(observations$COUNT) | observations$COUNT < 0,
           "Missing or invalid counts")
for (col in c("LW_A", "LW_B", "LENGTH_CONVERSION_FACTOR")) {
  values <- observations[[col]]
  check_rows(observations, !is.finite(values) | values <= 0,
             paste("Missing or invalid biomass coefficient", col))
}
weights <- with(observations, LW_A * (SIZE_ * LENGTH_CONVERSION_FACTOR)^LW_B)
check_rows(observations,
           !is.finite(weights) | !is.finite(weights * observations$COUNT),
           "Non-finite calculated biomass")
rm(survey_records, observations, weights)

for (reg in regions) {
  cat("\n------------------------------------------------------------\n")
  cat(sprintf("Processing Region: %s\n", reg))
  cat("------------------------------------------------------------\n")
  
  # Regional survey records
  reg_all <- nspc[REGION == reg]
  reg_obs <- nspc_obs[REGION == reg]
  
  if (nrow(reg_all) == 0) {
    stop("No records found for region: ", reg)
  }
  if (nrow(reg_obs) == 0) stop("No valid fish observations for region: ", reg)
  template_path <- sprintf("data/nSPC.site.abund.size.20092022.%s.rds", reg)
  if (!file.exists(template_path)) stop("Missing schema template: ", template_path)
  template <- readRDS(template_path)
  reference <- template[0, ]
  rm(template)
  
  # 1. Site Metadata: SITEVISITID, DEPTH, and site-level descriptors
  # Determine number of cylinders per site: n_distinct(REPLICATEID)
  site_cylinders <- reg_all[, .(n_cylinders = uniqueN(REPLICATEID)), by = SITEVISITID]
  
  # Site metadata with each unique DEPTH per site
  site_meta <- unique(reg_all[, .(
    METHOD = "nSPC",
    OBS_YEAR = as.integer(OBS_YEAR[1]),
    ANALYSIS_YEAR = as.character(ANALYSIS_YEAR[1]),
    ISLAND = as.character(ISLAND[1]),
    SEC_NAME = as.character(SEC_NAME[1]),
    REEF_ZONE = as.character(REEF_ZONE[1]),
    DEPTH_BIN = as.character(DEPTH_BIN[1]),
    LATITUDE = as.numeric(LATITUDE[1]),
    LONGITUDE = as.numeric(LONGITUDE[1]),
    DATE_ = DATE_[1]
  ), by = .(SITEVISITID, DEPTH)])
  
  site_meta <- merge(site_meta, site_cylinders, by = "SITEVISITID", all.x = TRUE)
  
  # 2. Species x Size_10cm grid for the region
  # Built from observed species-size combinations in the region with TAXONNAME
  sp_grid <- unique(reg_obs[, .(
    TAXONNAME = as.character(TAXONNAME[1])
  ), by = .(SPECIES, SIZE_10cm)])
  
  cat(sprintf("Region %s: %d sites, %d site-depth combinations, %d species-size grid rows\n",
              reg, uniqueN(site_meta$SITEVISITID), nrow(site_meta), nrow(sp_grid)))
  
  # 3. Calculate positive observations by SITEVISITID, SPECIES, SIZE_10cm
  cylinder_area <- pi * (7.5^2) # 176.7146 m^2
  
  # Calculate biomass per individual record
  reg_obs[, fish_wt := LW_A * ((SIZE_ * LENGTH_CONVERSION_FACTOR)^LW_B)]
  
  obs_agg <- reg_obs[, .(
    sum_count = sum(COUNT),
    sum_biom  = sum(fish_wt * COUNT)
  ), by = .(SITEVISITID, SPECIES, SIZE_10cm)]
  
  # 4. Expand full grid across site_meta and sp_grid
  # CJ (cross-join) site_meta and sp_grid
  cat("Expanding site-depth x species-size grid...\n")
  full_grid <- CJ(
    site_row = 1:nrow(site_meta),
    sp_row   = 1:nrow(sp_grid)
  )
  
  site_part <- site_meta[full_grid$site_row]
  sp_part   <- sp_grid[full_grid$sp_row]
  
  expanded <- cbind(site_part, sp_part)
  rm(full_grid, site_part, sp_part); gc()
  
  # 5. Merge observations into expanded grid
  cat("Merging observations and calculating abund.site and biom.site...\n")
  expanded <- merge(
    expanded,
    obs_agg,
    by = c("SITEVISITID", "SPECIES", "SIZE_10cm"),
    all.x = TRUE
  )
  
  expanded[is.na(sum_count), sum_count := 0]
  expanded[is.na(sum_biom), sum_biom := 0]
  
  # Density per depth entry = total / cylinder_area / n_cylinders
  expanded[, abund.site := (sum_count / cylinder_area) / n_cylinders]
  expanded[, biom.site  := (sum_biom / cylinder_area) / n_cylinders]
  
  # Factor conversions matching existing formats
  expanded[, `:=`(
    SITEVISITID   = as.integer(SITEVISITID),
    SPECIES       = as.character(SPECIES),
    SIZE_10cm     = as.character(SIZE_10cm),
    TAXONNAME     = preserve_factor(TAXONNAME, reference$TAXONNAME),
    METHOD        = preserve_factor(METHOD, reference$METHOD),
    OBS_YEAR      = as.integer(OBS_YEAR),
    ANALYSIS_YEAR = as.character(ANALYSIS_YEAR),
    ISLAND        = preserve_factor(ISLAND, reference$ISLAND),
    SEC_NAME      = as.character(SEC_NAME),
    REEF_ZONE     = preserve_factor(REEF_ZONE, reference$REEF_ZONE),
    DEPTH         = as.numeric(DEPTH),
    DEPTH_BIN     = preserve_factor(DEPTH_BIN, reference$DEPTH_BIN),
    LATITUDE      = as.numeric(LATITUDE),
    LONGITUDE     = as.numeric(LONGITUDE),
    DATE_         = as.POSIXct(DATE_)
  )]
  
  # Common 16 columns order
  cols_common <- c(
    "SITEVISITID", "SPECIES", "SIZE_10cm",
    "TAXONNAME", "METHOD", "OBS_YEAR", "ANALYSIS_YEAR", "ISLAND", "SEC_NAME",
    "REEF_ZONE", "DEPTH", "DEPTH_BIN", "LATITUDE", "LONGITUDE", "DATE_"
  )
  
  # 6. Save Abundance RDS
  abund_df <- expanded[, c(cols_common[1:3], "abund.site", cols_common[4:15]), with = FALSE]
  setDF(abund_df)
  if (!identical(names(abund_df), names(reference)) ||
      !identical(lapply(abund_df, class), lapply(reference, class))) {
    stop("Output schema does not match historical template for ", reg)
  }
  if (any(!is.finite(expanded$abund.site)) ||
      any(!is.finite(expanded$biom.site))) {
    stop("Non-finite output density for ", reg)
  }
  abund_out <- sprintf("data/nSPC.site.abund.size.20092024.%s.rds", reg)
  cat(sprintf("Saving %s (%d rows)...\n", abund_out, nrow(abund_df)))
  saveRDS(abund_df, file = abund_out)
  rm(abund_df); gc()
  
  # 7. Save Biomass RDS
  biom_df <- expanded[, c(cols_common[1:3], "biom.site", cols_common[4:15]), with = FALSE]
  setDF(biom_df)
  biom_out <- sprintf("data/nSPC.site.biom.size.20092024.%s.rds", reg)
  cat(sprintf("Saving %s (%d rows)...\n", biom_out, nrow(biom_df)))
  saveRDS(biom_df, file = biom_out)
  rm(biom_df, expanded); gc()
}

cat("\n=== Provisional nSPC Data Prep Finished; Normalization Still Unresolved ===\n")
