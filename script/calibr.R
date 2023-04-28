require(data.table) 
require(Calibr)
library(dplyr)
library(tidyr)
library(ggplot2)

rm(list = ls())

var = c("abund", "biom")[1]

belt <- readRDS(paste0("G:/ncrmp_spc/belt.site.", var, ".20002009.rds")) %>% 
  group_by(ISLAND, DEPTH_BIN, REEF_ZONE, METHOD, DATE_, LATITUDE, LONGITUDE, SPECIES) %>% 
  summarise(DENSITY = sum(!!sym(paste0(var, ".site"))), .groups = "drop") %>% 
  na.omit() %>% 
  mutate(PRESENCE = as.integer(DENSITY > 0),
         BLOCK = paste(ISLAND, DEPTH_BIN, REEF_ZONE, sep = "."),
         GROUP = SPECIES) %>% 
  select(DATE_, LATITUDE, LONGITUDE, BLOCK, GROUP, METHOD, DENSITY, PRESENCE)

spc = readRDS(paste0("data/nSPC.site.", var, ".size.20092022.MHI.rds")) %>% 
  group_by(ISLAND, DEPTH_BIN, REEF_ZONE, METHOD, DATE_, LATITUDE, LONGITUDE, SPECIES) %>% 
  summarise(DENSITY = sum(!!sym(paste0(var, ".site"))), .groups = "drop") %>% 
  na.omit() %>% 
  mutate(PRESENCE = as.integer(DENSITY > 0),
         BLOCK = paste(ISLAND, DEPTH_BIN, REEF_ZONE, sep = "."),
         GROUP = SPECIES) %>% 
  select(DATE_, LATITUDE, LONGITUDE, BLOCK, GROUP, METHOD, DENSITY, PRESENCE)

tow = readRDS(paste0("data/tow.segment.", var, ".size.20002017.MHI.rds")) %>% 
  subset(CENTROIDLON < 0) %>%
  mutate(DEPTH_BIN = case_when(
    DEPTH >= 0  & DEPTH <= 6 ~ "Shallow",
    DEPTH > 6  & DEPTH <= 18 ~ "Mid",
    DEPTH > 18 ~ "Deep",
    TRUE ~ ""),
    LONGITUDE = CENTROIDLON,
    LATITUDE = CENTROIDLAT) %>% 
  group_by(ISLAND, DEPTH_BIN, REEF_ZONE, METHOD, DATE_, LATITUDE, LONGITUDE, SPECIES) %>% 
  summarise(DENSITY = sum(!!sym(paste0(var, ".segment"))), .groups = "drop") %>% 
  na.omit() %>% 
  mutate(PRESENCE = as.integer(DENSITY > 0),
         BLOCK = paste(ISLAND, DEPTH_BIN, REEF_ZONE, sep = "."),
         GROUP = SPECIES)  %>% 
  select(DATE_, LATITUDE, LONGITUDE, BLOCK, GROUP, METHOD, DENSITY, PRESENCE)

calibr = c("spc_belt", "spc_tow")

for (c in 1:length(calibr)) {
  
  # c = 2
  
  if (calibr[c] == "spc_belt") set = rbind(spc, belt)
  if (calibr[c] == "spc_tow") set = rbind(spc, tow)

  set$REP <- sprintf("%04d", as.numeric(factor(paste(set$LATITUDE, set$LONGITUDE, set$DATE_, sep="_"))))
  set$REP = as.integer(set$REP)
  set = set[,c("BLOCK", "REP", "GROUP", "METHOD", "DENSITY", "PRESENCE")]
  set[] <- lapply(set, trimws)
  set$DENSITY = as.numeric(set$DENSITY)
  set$PRESENCE = as.numeric(set$PRESENCE)
  results <- run_calibr(set, std_method = "nSPC", stat_model = "GLMM", min_obs = 10)
  export_results(results, outdir = paste0("output/", calibr[c], "_", var))
  
}


