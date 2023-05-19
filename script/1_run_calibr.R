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

var = c("abund", "biom")
region = c("MHI", "MARIAN", "NWHI", "PRIAs", "SAMOA")
model = c("GLM", "GLMM")

spc_calibr = function(var, region, model){
  
  # var = "abund"
  # region = "MHI"
  # model = "GLM"
  
  belt <- readRDS(paste0("data/belt.site.", var, ".size.20002009.", region, ".rds")) %>% 
    group_by(SITEVISITID, SPECIES, METHOD, OBS_YEAR, ISLAND, REEF_ZONE, DEPTH_BIN, LATITUDE, LONGITUDE, DATE_, n.transect) %>% # aggregate size bins
    summarise(DENSITY = sum(!!sym(paste0(var, ".site")))) %>% 
    mutate(PRESENCE = as.integer(DENSITY > 0)) %>% 
    group_by(SITEVISITID) %>%
    mutate(unique_transect_count = n_distinct(n.transect)) %>%
    ungroup() %>% 
    mutate(DENSITY = DENSITY/unique_transect_count,
           PRESENCE = PRESENCE/unique_transect_count, 
           BLOCK = paste(ISLAND, DEPTH_BIN, REEF_ZONE, sep = "."),
           # BLOCK = paste(OBS_YEAR, ISLAND, DEPTH_BIN, REEF_ZONE, sep = "."),
           GROUP = SPECIES) %>% 
    select(DATE_, LATITUDE, LONGITUDE, BLOCK, GROUP, METHOD, DENSITY, PRESENCE)
  
  spc = readRDS(paste0("data/nSPC.site.", var, ".size.20092022.", region, ".rds")) %>% 
    group_by(SITEVISITID, SPECIES, METHOD, OBS_YEAR, ISLAND, REEF_ZONE, DEPTH_BIN, LATITUDE, LONGITUDE, DATE_) %>% # aggregate size bins
    summarise(DENSITY = sum(!!sym(paste0(var, ".site")))) %>% 
    ungroup() %>%
    mutate(PRESENCE = as.integer(DENSITY > 0), 
           BLOCK = paste(ISLAND, DEPTH_BIN, REEF_ZONE, sep = "."),
           # BLOCK = paste(OBS_YEAR, ISLAND, DEPTH_BIN, REEF_ZONE, sep = "."),
           GROUP = SPECIES) %>% 
    select(DATE_, LATITUDE, LONGITUDE, BLOCK, GROUP, METHOD, DENSITY, PRESENCE)
  
  tow = readRDS(paste0("data/tow.segment.", var, ".size.20002017.", region, ".rds")) %>% 
    subset(CENTROIDLON != 0 & SIZE_10cm != "(40,50]") %>% # remove size bins < 50 cm
    mutate(DEPTH_BIN = case_when(
      DEPTH >= 0  & DEPTH <= 6 ~ "Shallow",
      DEPTH > 6  & DEPTH <= 18 ~ "Mid",
      DEPTH > 18 ~ "Deep",
      TRUE ~ ""),
      LONGITUDE = CENTROIDLON,
      LATITUDE = CENTROIDLAT) %>% 
    group_by(TOWID, SEGMENTID, SPECIES, METHOD, OBS_YEAR, ISLAND, REEF_ZONE, DEPTH_BIN, LATITUDE, LONGITUDE, DATE_) %>% # aggregate size bins
    summarise(DENSITY = sum(!!sym(paste0(var, ".segment")))) %>% 
    mutate(PRESENCE = as.integer(DENSITY > 0)) %>% 
    ungroup() %>% 
    group_by(TOWID) %>%
    mutate(unique_segment_count = n_distinct(SEGMENTID)) %>% 
    ungroup() %>% 
    group_by(TOWID, SPECIES, METHOD, OBS_YEAR, ISLAND, REEF_ZONE, DEPTH_BIN, LATITUDE, LONGITUDE, DATE_) %>% # average by unique_segment_count
    summarise(DENSITY = DENSITY/unique_segment_count,
              PRESENCE = PRESENCE/unique_segment_count) %>% 
    mutate(
      BLOCK = paste(ISLAND, DEPTH_BIN, REEF_ZONE, sep = "."),
      # BLOCK = paste(OBS_YEAR, ISLAND, DEPTH_BIN, REEF_ZONE, sep = "."),
      GROUP = SPECIES) %>% 
    ungroup() %>% 
    select(DATE_, LATITUDE, LONGITUDE, BLOCK, GROUP, METHOD, DENSITY, PRESENCE)
  
  calibr = c("spc_belt", "spc_tow")
  
  for (c in 1:length(calibr)) {
    
    # c = 1
    
    if (calibr[c] == "spc_belt") set = rbind(spc, belt)
    if (calibr[c] == "spc_tow") set = rbind(spc, tow)
    
    set$REP <- sprintf("%04d", as.numeric(factor(paste(set$LATITUDE, set$LONGITUDE, set$DATE_, sep="_"))))
    set$REP = as.integer(set$REP)
    set = set[,c("BLOCK", "REP", "GROUP", "METHOD", "DENSITY", "PRESENCE")]
    set[] <- lapply(set, trimws)
    set$DENSITY = as.numeric(set$DENSITY)
    set$PRESENCE = as.numeric(set$PRESENCE)
    results <- run_calibr(set, std_method = "nSPC", stat_model = model)
    if(!file.exists(paste0("output/", calibr[c], "_", var, "_", region, "_", model))) dir.create(paste0("output/", calibr[c], "_", var, "_", region, "_", model))
    export_results(results, outdir = paste0("output/", calibr[c], "_", var, "_", region, "_", model))
    
  }
  
  
}

spc_calibr("abund", "MHI", "GLM")
spc_calibr("abund", "MHI", "GLMM")
spc_calibr("biom", "MHI", "GLM")
spc_calibr("biom", "MHI", "GLMM")

spc_calibr("abund", "SAMOA", "GLM")
spc_calibr("abund", "SAMOA", "GLMM")
spc_calibr("biom", "SAMOA", "GLM")
spc_calibr("biom", "SAMOA", "GLMM")

spc_calibr("abund", "MARIAN", "GLM")
spc_calibr("abund", "MARIAN", "GLMM")
spc_calibr("biom", "MARIAN", "GLM")
spc_calibr("biom", "MARIAN", "GLMM")

spc_calibr("abund", "NWHI", "GLM")
spc_calibr("abund", "NWHI", "GLMM")
spc_calibr("biom", "NWHI", "GLM")
spc_calibr("biom", "NWHI", "GLMM")

spc_calibr("abund", "PRIAs", "GLM")
spc_calibr("abund", "PRIAs", "GLMM")
spc_calibr("biom", "PRIAs", "GLM")
spc_calibr("biom", "PRIAs", "GLMM")