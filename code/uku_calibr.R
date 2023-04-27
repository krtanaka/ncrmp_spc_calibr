require(data.table) 
require(Calibr)
library(dplyr)
library(tidyr)
library(ggplot2)

rm(list = ls())

belt = readRDS("data/belt.site.abund.size.20002009.MHI.rds") %>% 
  group_by(ISLAND, DEPTH_BIN, REEF_ZONE, METHOD, DATE_, LATITUDE, LONGITUDE, SPECIES) %>% 
  summarise(DENSITY = sum(abund.site)) %>% 
  na.omit() %>% 
  mutate(PRESENCE = ifelse(DENSITY > 0, 1, 0),
         BLOCK = paste(ISLAND, DEPTH_BIN, REEF_ZONE, sep = "."),
         GROUP = SPECIES) %>% 
  as.data.frame()

belt = belt[,c("DATE_", "LATITUDE", "LONGITUDE", "BLOCK", "GROUP", "METHOD", "DENSITY", "PRESENCE")]

spc = readRDS("data/nSPC.site.abund.size.20092022.MHI.rds") %>% 
  group_by(ISLAND, DEPTH_BIN, REEF_ZONE, METHOD, DATE_, LATITUDE, LONGITUDE, SPECIES) %>% 
  summarise(DENSITY = sum(abund.site)) %>% 
  na.omit() %>% 
  mutate(PRESENCE = ifelse(DENSITY > 0, 1, 0),
         BLOCK = paste(ISLAND, DEPTH_BIN, REEF_ZONE, sep = "."),
         GROUP = SPECIES) %>% 
  as.data.frame()

spc = spc[,c("DATE_", "LATITUDE", "LONGITUDE", "BLOCK", "GROUP", "METHOD", "DENSITY", "PRESENCE")]

tow = readRDS("data/tow.segment.abund.size.20002017.MHI.rds") %>% 
  mutate(
    DEPTH_BIN = case_when(
      DEPTH >= 0  & DEPTH <= 6 ~ "Shallow",
      DEPTH > 6  & DEPTH <= 18 ~ "Mid",
      DEPTH > 18 ~ "Deep",
      TRUE ~ ""),
    LONGITUDE = CENTROIDLON,
    LATITUDE = CENTROIDLAT) %>% 
  group_by(ISLAND, DEPTH_BIN, REEF_ZONE, METHOD, DATE_, LATITUDE, LONGITUDE, SPECIES) %>% 
  summarise(DENSITY = sum(abund.segment)) %>%
  na.omit() %>% 
  mutate(PRESENCE = ifelse(DENSITY > 0, 1, 0),
         BLOCK = paste(ISLAND, DEPTH_BIN, REEF_ZONE, sep = "."),
         GROUP = SPECIES) %>% 
  as.data.frame()

tow = tow[,c("DATE_", "LATITUDE", "LONGITUDE", "BLOCK", "GROUP", "METHOD", "DENSITY", "PRESENCE")]

set = rbind(spc, belt)
set$REP <- sprintf("%04d", as.numeric(factor(paste(set$LATITUDE, set$LONGITUDE, set$DATE_, sep="_"))))
set$REP = as.integer(set$REP)
set = set[,c("BLOCK", "REP", "GROUP", "METHOD", "DENSITY", "PRESENCE")]
set = set %>% as.data.frame() %>% na.omit() %>% subset(GROUP == "APVI")
# set$DENSITY = round(set$DENSITY*1000, 1)


head(SMALL_UNPAIR)
head(set)

results <- run_calibr(SMALL_UNPAIR, std_method="nSPC", stat_model="GLM", min_obs = 10)
export_results(results)

results$LGROUP