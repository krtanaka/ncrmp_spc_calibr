library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(readr)
library(tidyverse)

rm(list = ls())

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

folders <- list.dirs(path = "output", full.names = TRUE, recursive = FALSE)

read_summary_table <- function(folder_path) {
  file_path <- file.path(folder_path, "summary_table.csv")
  if (file.exists(file_path)) {
    df <- read.csv(file_path)
    return(df)
  } else {
    message(paste("File not found:", file_path))
    return(NULL)
  }
}

df_list <- lapply(folders, read_summary_table)
names(df_list) <- basename(folders)

species = "APVI"

png("output/calibr_APVI_coef.png", units = "in", height = 3, width = 7, res = 500)

(df <- bind_rows(df_list, .id = "folder")  %>% 
    separate(folder, into = c("spc", "belt_tow", "var", "region", "model"), sep = "_") %>%
    rename_all(tolower) %>% 
    filter(group == species, model == "GLM", method != "1_nSPC") %>%
    ggplot(aes(region, gcf.pos, color = var)) +
    geom_point(position = position_dodge(width = 0.3), size = 2) +
    geom_errorbar(aes(ymin = gcf.pos_2.5, ymax = gcf.pos_95), 
                  position = position_dodge(width = 0.3), width = 0, show.legend = F) +
    facet_wrap(~belt_tow, scales = "free") + 
    scale_color_discrete(species) + 
    labs(x = "Region", y = "Gear Calibration Factor"))

dev.off()
