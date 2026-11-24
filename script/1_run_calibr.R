# install.packages("devtools")
# library(devtools)
# install_github("marcnadon/Calibr")

# # Marc's example
# SET <- SMALL_UNPAIR
# results <- run_calibr(SET,std_method="nSPC",stat_model="GLM")
# export_results(results)

require(data.table) 
require(Calibr)
library(dplyr)
library(tidyr)
library(ggplot2)

rm(list = ls())

source("script/nspc_summary.R")

calibr_args <- commandArgs(trailingOnly = TRUE)
if (length(calibr_args) > 5L) {
  stop("Usage: 1_run_calibr.R [GROUP [METRIC [REGION [MODEL [METHOD]]]]]",
       call. = FALSE)
}
calibr_groups <- if (length(calibr_args)) calibr_args[[1]] else NULL
metric_filter <- if (length(calibr_args) >= 2L) calibr_args[[2]] else NULL
region_filter <- if (length(calibr_args) >= 3L) calibr_args[[3]] else NULL
model_filter <- if (length(calibr_args) >= 4L) calibr_args[[4]] else NULL
method_filter <- if (length(calibr_args) >= 5L) calibr_args[[5]] else NULL
output_root <- Sys.getenv("CALIBR_OUTPUT_DIR", unset = "output")
if (!dir.exists(output_root)) {
  stop("Calibration output directory does not exist: ", output_root, call. = FALSE)
}
if (!is.null(metric_filter) && !metric_filter %in% c("abund", "biom")) {
  stop("Unsupported metric: ", metric_filter, call. = FALSE)
}
if (!is.null(region_filter) &&
    !region_filter %in% c("MHI", "MARIAN", "NWHI", "PRIAs", "SAMOA")) {
  stop("Unsupported region: ", region_filter, call. = FALSE)
}
if (!is.null(model_filter) && !model_filter %in% c("GLM", "GLMM")) {
  stop("Unsupported model: ", model_filter, call. = FALSE)
}
if (!is.null(method_filter) && !method_filter %in% c("belt", "tow")) {
  stop("Unsupported method: ", method_filter, call. = FALSE)
}
if (is.null(model_filter) || model_filter == "GLMM") {
  glmm_source <- paste(deparse(body(getFromNamespace("gcf_glmm", "Calibr"))),
                       collapse = " ")
  if (grepl('gsub(paste0(":1_",', glmm_source, fixed = TRUE)) {
    stop("Installed Calibr has a GLMM group-name bug; see README for the fix.",
         call. = FALSE)
  }
}

merge_calibr_exports <- function(results, outdir, groups) {
  temp_dir <- tempfile("calibr-export-")
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)
  export_results(results, outdir = temp_dir)

  for (filename in c("summary_table.csv", "REP_summary.csv")) {
    old_path <- file.path(outdir, filename)
    new_path <- file.path(temp_dir, filename)
    if (!file.exists(old_path) || !file.exists(new_path)) {
      stop("Missing calibration export: ", old_path, " or ", new_path,
           call. = FALSE)
    }

    old_lines <- readLines(old_path)
    new_lines <- readLines(new_path)
    if (!identical(old_lines[1], new_lines[1])) {
      stop("Mismatched headers in calibration export: ", filename,
           call. = FALSE)
    }

    for (group in groups) {
      old_rows <- which(startsWith(old_lines, paste0('"', group, '",')))
      new_rows <- which(startsWith(new_lines, paste0('"', group, '",')))
      if (!length(old_rows) || !length(new_rows)) {
        stop("Missing calibration rows for ", group, " in ", filename,
             call. = FALSE)
      }
      if (length(old_rows) == length(new_rows)) {
        old_lines[old_rows] <- new_lines[new_rows]
      } else {
        first_old_row <- min(old_rows)
        insert_at <- if (first_old_row > 1L) {
          sum(!startsWith(old_lines[seq_len(first_old_row - 1L)],
                          paste0('"', group, '",')))
        } else {
          0L
        }
        old_lines <- old_lines[-old_rows]
        old_lines <- append(old_lines, new_lines[new_rows], after = insert_at)
      }
    }
    writeLines(old_lines, old_path, useBytes = TRUE)
  }

  old_path <- file.path(outdir, "reeffish_grouped.Rds")
  new_path <- file.path(temp_dir, "reeffish_grouped.Rds")
  if (!file.exists(old_path) || !file.exists(new_path)) {
    stop("Missing grouped calibration export: ", old_path, " or ", new_path,
         call. = FALSE)
  }
  old_results <- readRDS(old_path)
  new_results <- readRDS(new_path)
  for (group in groups) {
    if (is.null(new_results[[group]])) {
      stop("Missing grouped calibration results for ", group, call. = FALSE)
    }
    old_results[[group]] <- new_results[[group]]
  }
  saveRDS(old_results, old_path)
}

spc_calibr = function(var, region, model, groups = calibr_groups) {
  if (!var %in% c("abund", "biom")) stop("Unsupported metric: ", var)
  if (!region %in% c("MHI", "MARIAN", "NWHI", "PRIAs", "SAMOA")) {
    stop("Unsupported region: ", region)
  }
  if (!model %in% c("GLM", "GLMM")) stop("Unsupported model: ", model)
  
  # belt transect data
  belt <- readRDS(paste0("data/belt.site.", var, ".size.20002009.", region, ".rds")) %>%
    filter(SIZE_10cm != "(0,10]" & SIZE_10cm != "(10,20]") %>%
    group_by(SITEVISITID, SPECIES, METHOD, OBS_YEAR, ISLAND, REEF_ZONE, DEPTH_BIN, LATITUDE, LONGITUDE, DATE_, n.transect) %>% 
    summarise(DENSITY = sum(!!sym(paste0(var, ".site")))) %>% 
    mutate(PRESENCE = as.integer(DENSITY > 0),
           BLOCK = paste(ISLAND, DEPTH_BIN, REEF_ZONE, sep = "."),
           # BLOCK = paste(OBS_YEAR, ISLAND, DEPTH_BIN, REEF_ZONE, sep = "."),
           GROUP = SPECIES
    ) %>% 
    mutate(
      DENSITY = DENSITY/n.transect,
      PRESENCE = PRESENCE/n.transect,
    ) %>%
    ungroup() %>%
    select(DATE_, LATITUDE, LONGITUDE, BLOCK, GROUP, METHOD, DENSITY, PRESENCE)
  
  # towed diver data
  tow = readRDS(paste0("data/tow.segment.", var, ".size.20002017.", region, ".rds")) %>% 
    subset(CENTROIDLON != 0 & SIZE_10cm != "(40,50]") %>%
    mutate(DEPTH_BIN = case_when(
      DEPTH >= 0  & DEPTH <= 6 ~ "Shallow",
      DEPTH > 6  & DEPTH <= 18 ~ "Mid",
      DEPTH > 18 ~ "Deep",
      TRUE ~ ""),
      LONGITUDE = CENTROIDLON,
      LATITUDE = CENTROIDLAT) %>% 
    group_by(TOWID, SEGMENTID, SPECIES, METHOD, OBS_YEAR, ISLAND, REEF_ZONE, DEPTH_BIN, LATITUDE, LONGITUDE, DATE_) %>% 
    summarise(DENSITY = sum(!!sym(paste0(var, ".segment")))) %>% 
    mutate(
      PRESENCE = as.integer(DENSITY > 0),
      BLOCK = paste(ISLAND, DEPTH_BIN, REEF_ZONE, sep = "."),
      # BLOCK = paste(OBS_YEAR, ISLAND, DEPTH_BIN, REEF_ZONE, sep = "."),
      GROUP = SPECIES
    ) %>% 
    ungroup() %>% 
    group_by(TOWID) %>%
    mutate(
      unique_segment_count = n_distinct(SEGMENTID),
      DENSITY = DENSITY/unique_segment_count,
      PRESENCE = PRESENCE/unique_segment_count) %>%
    ungroup() %>% 
    select(DATE_, LATITUDE, LONGITUDE, BLOCK, GROUP, METHOD, DENSITY, PRESENCE)
  
  # Use the I-only summary because it is the calibration reference dataset.
  spc_density <- if (var == "abund") "DENSITY_M2" else "BIOMASS_G_M2"
  spc <- read_nspc_summary(
    "I", region,
    columns = c("METHOD", "ISLAND", "REEF_ZONE", "DEPTH_BIN", "DATE_",
                "LATITUDE", "LONGITUDE", "FISH_GROUP", "PRESENT",
                spc_density)
  ) %>%
    transmute(
      DATE_, LATITUDE, LONGITUDE,
      BLOCK = paste(ISLAND, DEPTH_BIN, REEF_ZONE, sep = "."),
      GROUP = FISH_GROUP, METHOD,
      DENSITY = .data[[spc_density]],
      PRESENCE = as.numeric(PRESENT)
    )
  
  calibr = c("spc_belt", "spc_tow")
  if (!is.null(method_filter)) {
    calibr <- paste0("spc_", method_filter)
  }
  
  for (c in 1:length(calibr)) {
    
    # c = 2
    
    if (calibr[c] == "spc_belt") set = rbind(spc, belt)
    if (calibr[c] == "spc_tow") set = rbind(spc, tow)

    if (!is.null(groups)) {
      set <- set %>% filter(GROUP %in% groups)
      if (!nrow(set)) {
        stop("No calibration records found for: ",
             paste(groups, collapse = ", "), call. = FALSE)
      }
    }
    
    set$REP <- sprintf("%04d", as.numeric(factor(paste(set$LATITUDE, set$LONGITUDE, set$DATE_, sep = "_"))))
    set$REP = as.integer(set$REP)
    set = set[,c("BLOCK", "REP", "GROUP", "METHOD", "DENSITY", "PRESENCE")]
    set[] <- lapply(set, trimws)
    set$DENSITY = as.numeric(set$DENSITY)
    set$PRESENCE = as.numeric(set$PRESENCE)
    
    results <- run_calibr(set, std_method = "nSPC", stat_model = model)
    
    outdir <- file.path(output_root, paste0(calibr[c], "_", var, "_", region, "_", model))
    if (!file.exists(outdir)) dir.create(outdir)
    if (is.null(groups)) {
      export_results(results, outdir = outdir)
    } else {
      merge_calibr_exports(results, outdir, groups)
    }
    
  }

}

run_calibration <- function(var, region, model) {
  if (!is.null(metric_filter) && var != metric_filter) return(invisible(NULL))
  if (!is.null(region_filter) && region != region_filter) return(invisible(NULL))
  if (!is.null(model_filter) && model != model_filter) return(invisible(NULL))
  spc_calibr(var, region, model)
}

run_calibration("abund", "MHI", "GLM")
run_calibration("abund", "MHI", "GLMM")
run_calibration("biom", "MHI", "GLM")
run_calibration("biom", "MHI", "GLMM")

run_calibration("abund", "SAMOA", "GLM")
run_calibration("abund", "SAMOA", "GLMM")
run_calibration("biom", "SAMOA", "GLM")
run_calibration("biom", "SAMOA", "GLMM")

run_calibration("abund", "MARIAN", "GLM")
run_calibration("abund", "MARIAN", "GLMM")
run_calibration("biom", "MARIAN", "GLM")
run_calibration("biom", "MARIAN", "GLMM")

run_calibration("abund", "NWHI", "GLM")
run_calibration("abund", "NWHI", "GLMM")
run_calibration("biom", "NWHI", "GLM")
run_calibration("biom", "NWHI", "GLMM")

run_calibration("abund", "PRIAs", "GLM")
run_calibration("abund", "PRIAs", "GLMM")
run_calibration("biom", "PRIAs", "GLM")
run_calibration("biom", "PRIAs", "GLMM")