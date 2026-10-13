# ==============================================================================
# script/0_prep_nspc_data.R
#
# Process raw survey data from fish-paste (through 2024) and generate regional
# nSPC site-level abundance and biomass RDS files preserving exact schema,
# column types, and factor levels for downstream calibration models.
# ==============================================================================

library(data.table)

cat("=== Starting nSPC Data Prep (2009-2024) ===\n")

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

cat("Loading ALL_REA_FISH_RAW...\n")
load(raw_data_path) # loads `df`
setDT(df)

cat("Filtering valid nSPC survey records (2009-2024)...\n")
# Filter nSPC surveys
nspc <- df[
  METHOD %in% c("nSPC", "nSPC-CCR") &
  OBS_YEAR > 2008 &
  (is.na(EXCLUDE_FLAG) | EXCLUDE_FLAG == 0) &
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

# Standardize species names: convert CHPA to CHAG
nspc[SPECIES == "CHPA", SPECIES := "CHAG"]

# Calculate 10-cm size bins: (0,10], (10,20], ..., (290,300]
nspc[, SIZE_10cm := as.character(cut(SIZE_, breaks = seq(0, 300, 10)))]

# Filter valid fish observation types for abundance/biomass
# (U: unconstrained, I: instantaneous, N: non-instantaneous)
nspc_obs <- nspc[OBS_TYPE %in% c("U", "I", "N") & !is.na(SPECIES) & !is.na(SIZE_10cm)]

# Target regions
regions <- c("MARIAN", "MHI", "NWHI", "PRIAs", "SAMOA")

for (reg in regions) {
  cat("\n------------------------------------------------------------\n")
  cat(sprintf("Processing Region: %s\n", reg))
  cat("------------------------------------------------------------\n")
  
  # Regional survey records
  reg_all <- nspc[REGION == reg]
  reg_obs <- nspc_obs[REGION == reg]
  
  if (nrow(reg_all) == 0) {
    warning("No records found for region: ", reg)
    next
  }
  
  # 1. Site Metadata: SITEVISITID, DEPTH, and site-level descriptors
  # Determine number of cylinders per site: n_distinct(REPLICATEID)
  site_cylinders <- reg_all[, .(n_cylinders = uniqueN(REPLICATEID)), by = SITEVISITID]
  
  # Site metadata with each unique DEPTH per site
  site_meta <- unique(reg_all[, .(
    SITEVISITID,
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
    SPECIES = as.character(SPECIES),
    SIZE_10cm = as.character(SIZE_10cm),
    TAXONNAME = as.character(TAXONNAME[1])
  ), by = .(SPECIES, SIZE_10cm)])
  
  cat(sprintf("Region %s: %d sites, %d site-depth combinations, %d species-size grid rows\n",
              reg, uniqueN(site_meta$SITEVISITID), nrow(site_meta), nrow(sp_grid)))
  
  # 3. Calculate positive observations by SITEVISITID, SPECIES, SIZE_10cm
  cylinder_area <- pi * (7.5^2) # 176.7146 m^2
  
  # Calculate biomass per individual record
  reg_obs[, fish_wt := LW_A * ((SIZE_ * LENGTH_CONVERSION_FACTOR)^LW_B)]
  
  obs_agg <- reg_obs[, .(
    sum_count = sum(COUNT, na.rm = TRUE),
    sum_biom  = sum(fish_wt * COUNT, na.rm = TRUE)
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
    TAXONNAME     = as.factor(TAXONNAME),
    METHOD        = as.factor(METHOD),
    OBS_YEAR      = as.integer(OBS_YEAR),
    ANALYSIS_YEAR = as.character(ANALYSIS_YEAR),
    ISLAND        = as.factor(ISLAND),
    SEC_NAME      = as.character(SEC_NAME),
    REEF_ZONE     = as.factor(REEF_ZONE),
    DEPTH         = as.numeric(DEPTH),
    DEPTH_BIN     = as.factor(DEPTH_BIN),
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

cat("\n=== nSPC Data Prep Finished Successfully ===\n")
