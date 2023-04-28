library(readr)

rm(list = ls())

var = c("abund", "biom")[1]

# breakdown of replicate sample sizes
# NREP_TOTAL: Total number of replicates in the data set (both methods combined).This excludes replicates from BLOCKs where a specific GROUP was not found.
# NREP_STD_METHOD: Total number of replicates in the data set for the standard method.
# NREP_SEC_METHOD: Total number of replicates in the data set for the secondary method.
# POSREP_TOTAL: Total number of replicates where a specific GROUP was observed (positive, non-zero count).
# POSREP_STD_METHOD: Total number of replicates from standard method where a specific GROUP was observed (positive, non-zero count).
# POSREP_SEC_METHOD: Total number of replicates from secondary method where a specific GROUP was observed (positive, non-zero count).

read_csv(paste0("output/spc_belt_", var, "/REP_summary.csv")) %>% subset(GROUP == "APVI")
read_csv(paste0("output/spc_tow_", var, "/REP_summary.csv")) %>% subset(GROUP == "APVI")

# final results
# METHOD: Sampling method used.
# GCF.PRES: Gear calibration factor for presence-absence data.
# Convert probability of observation of secondary method following this equation: 
# Prob_M1=inv.logit(logit(Prob_M2-GCF.PRES)))
# GCF.PRES_2.5: Lower bound of 95% probability interval of GCF.PRES.
# GCF.PRES_95: Upper bound of 95% probability interval of GCF.PRES.
# GCF.POS: Gear calibration factor for positive-only data. Convert abundance metric of positive-only data for secondary method following this equation: 
# Abund_M1=Abund_M2/GCF.POS
# GCF.POS_2.5: Lower bound of 95% probability interval of GCF.PRES.
# GCF.POS_95: Upper bound of 95% probability interval of GCF.PRES.
# PRES: Probability of observing a specific species for each method in the provided dataset.
# PRES.CAL: Calibrated probability of observing a species by method (used as a check on the success of the standardization procedure).
# POS: Abundance of a specific species for each method.
# POS.CAL: Calibrated abundance of a specific species for each method (used as a check on the success of the standardization procedure).
# OPUE: Obervation Per Unit Effort obtained by multiplying prob. of sighting and abundance (PRES x POS) in provided data set.
# OPUE.CAL: Obervation Per Unit Effort obtained by multiplying calibrated prob. of sighting and abundance (PRES.CAL x POS.CAL) in provided data set.

calibr_belt = read_csv(paste0("output/spc_belt_", var, "/summary_table.csv")) %>% subset(GROUP == "APVI"); calibr_belt = calibr_belt$POS.GCF
calibr_tow = read_csv(paste0("output/spc_tow_", var, "/summary_table.csv")) %>% subset(GROUP == "APVI"); calibr_tow = calibr_tow$POS.GCF

belt <- readRDS(paste0("data/belt.site.", var, ".size.20002009.MHI.rds")) %>% 
  subset(SPECIES == "APVI") %>% 
  group_by(ISLAND, DEPTH, METHOD, DATE_, LATITUDE, LONGITUDE) %>% 
  summarise(DENSITY = sum(!!sym(paste0(var, ".site"))), .groups = "drop") %>% 
  mutate(DENSITY = (DENSITY / calibr_belt[1]))

spc = readRDS(paste0("data/nSPC.site.", var, ".size.20092022.MHI.rds")) %>% 
  subset(SPECIES == "APVI") %>% 
  group_by(ISLAND, DEPTH, METHOD, DATE_, LATITUDE, LONGITUDE) %>% 
  summarise(DENSITY = sum(!!sym(paste0(var, ".site"))), .groups = "drop") %>% 
  mutate(DENSITY = DENSITY)

tow = readRDS(paste0("data/tow.segment.", var, ".size.20002017.MHI.rds")) %>% 
  mutate(LONGITUDE = CENTROIDLON,
         LATITUDE = CENTROIDLAT) %>% 
  group_by(ISLAND, DEPTH, METHOD, DATE_, LATITUDE, LONGITUDE) %>% 
  summarise(DENSITY = sum(!!sym(paste0(var, ".segment"))), .groups = "drop") %>% 
  mutate(DENSITY = (DENSITY / calibr_tow[1]))

df = rbind(spc, belt, tow)

df %>% 
  ggplot(aes(LONGITUDE, LATITUDE, size = DENSITY)) + 
  geom_point()
