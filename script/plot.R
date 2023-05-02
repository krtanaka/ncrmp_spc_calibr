library(readr)
library(lubridate)
library(colorRamps)
library(ggplot2)
library(dplyr)

rm(list = ls())

region = c("MHI", "MARIAN", "NWHI", "PRIAs", "SAMOA")
var = c("abund", "biom")[2]
species = "APVI"

df = NULL

for (r in 1:length(region)) {
  
  calibr_belt = read_csv(paste0("output/spc_belt_", var, "_", region[r], "_GLMM/summary_table.csv")) %>% 
    filter(GROUP == species) %>% 
    filter(METHOD != "1_nSPC")
  
  calibr_tow = read_csv(paste0("output/spc_tow_", var, "_", region[r], "_GLMM/summary_table.csv")) %>% 
    filter(GROUP == species) %>% 
    filter(METHOD != "1_nSPC")
  
  calibr_tow = calibr_tow$POS.GCF
  calibr_belt = calibr_belt$POS.GCF
  
  belt <- readRDS(paste0("data/belt.site.", var, ".size.20002009.", region[r], ".rds")) %>% 
    group_by(ISLAND, DEPTH, METHOD, DATE_, LATITUDE, LONGITUDE, SPECIES) %>% 
    summarise(DENSITY = sum(!!sym(paste0(var, ".site"))), .groups = "drop")  %>% 
    filter(SPECIES == species) %>%
    mutate(DENSITY = (DENSITY / calibr_belt))
  
  spc = readRDS(paste0("data/nSPC.site.", var, ".size.20092022.", region[r], ".rds"))  %>% 
    group_by(ISLAND, DEPTH, METHOD, DATE_, LATITUDE, LONGITUDE, SPECIES) %>% 
    summarise(DENSITY = sum(!!sym(paste0(var, ".site"))), .groups = "drop") %>% 
    filter(SPECIES == species)
  
  tow = readRDS(paste0("data/tow.segment.", var, ".size.20002017.", region[r], ".rds"))  %>% 
    filter(CENTROIDLON != 0) %>%
    mutate(LONGITUDE = CENTROIDLON,
           LATITUDE = CENTROIDLAT) %>% 
    group_by(ISLAND, DEPTH, METHOD, DATE_, LATITUDE, LONGITUDE, SPECIES) %>% 
    summarise(DENSITY = sum(!!sym(paste0(var, ".segment"))), .groups = "drop") %>% 
    filter(SPECIES == species) %>%
    mutate(DENSITY = (DENSITY / calibr_tow))
  
  df_i = rbind(spc, belt, tow)
  df_i$region = region[r]
  
  df = rbind(df, df_i)
  
}

df <- df %>%
  rename_all(tolower) %>% 
  mutate(year = year(date_),
         month = month(date_),
         day = day(date_),
         longitude = ifelse(longitude < 0, longitude + 360, longitude),
         density = density*100)

save(df, file = paste0("output/calibr_", species, "_", var, ".RData"))

df %>% 
  mutate(longitude = round(longitude, 1), 
         latitude = round(latitude, 1)) %>% 
  group_by(method, longitude, latitude, region) %>%
  summarise(density = mean(density)) %>%
  ggplot(aes(longitude, latitude)) + 
  geom_point(aes(size = density, fill = density, color = density), shape = 21, alpha = 0.7) +
  scale_color_gradientn(colours = matlab.like(100), guide = "legend", trans = "sqrt") +
  scale_fill_gradientn(colours = matlab.like(100), guide = "legend", trans = "sqrt") +
  facet_wrap(method ~ region, scales = "free") +
  labs(x = expression(paste("Longitude ", degree, "W", sep = "")),
       y = expression(paste("Latitude ", degree, "N", sep = ""))) +
  guides(color = guide_legend(expression("Individuals per 100" ~ m^2~"")), 
         fill = guide_legend(expression("Individuals per 100" ~ m^2~"")),
         size = guide_legend(expression("Individuals per 100" ~ m^2~"")))

df %>% 
  filter(region == "MHI") %>% 
  mutate(density = density*100,
         longitude = round(longitude, 1), 
         latitude = round(latitude, 1)) %>% 
  group_by(method, longitude, latitude, year) %>%
  summarise(density = mean(density)) %>%
  ggplot(aes(longitude, latitude)) + 
  geom_point(aes(size = density, fill = density, color = density), shape = 21, alpha = 0.7) +
  scale_color_gradientn(colours = matlab.like(100), guide = "legend", trans = "sqrt") +
  scale_fill_gradientn(colours = matlab.like(100), guide = "legend", trans = "sqrt") +
  facet_grid(method ~ year) +
  labs(x = expression(paste("Longitude ", degree, "W", sep = "")),
       y = expression(paste("Latitude ", degree, "N", sep = ""))) +
  guides(color = guide_legend(expression("Individuals per 100" ~ m^2~"")), 
         fill = guide_legend(expression("Individuals per 100" ~ m^2~"")),
         size = guide_legend(expression("Individuals per 100" ~ m^2~"")))

df %>% 
  mutate(density = density*100) %>% 
  ggplot(aes(depth, density)) + 
  geom_point(aes(size = density, fill = density, color = density), shape = 21, alpha = 0.8) +
  scale_color_gradientn(colours = matlab.like(100), guide = "legend", trans = "sqrt") +
  scale_fill_gradientn(colours = matlab.like(100), guide = "legend", trans = "sqrt") +
  facet_grid( ~ method) +
  guides(color = guide_legend(expression("Individuals per 100" ~ m^2~"")), 
         fill = guide_legend(expression("Individuals per 100" ~ m^2~"")),
         size = guide_legend(expression("Individuals per 100" ~ m^2~"")))

df %>%
  mutate(YEAR = format(date_, "%Y"),
         density = density*100) %>% 
  group_by(year, method, region) %>%
  summarise(mean_density = mean(density), se_density = sd(density)/sqrt(n())) %>%
  mutate(mean_density = ifelse(mean_density == 0, NA, mean_density),
         se_density = ifelse(se_density == 0, NA, se_density)) %>% 
  ggplot(aes(x = year, y = mean_density, color = method, group = method)) +
  geom_errorbar(aes(ymin = mean_density - se_density, ymax = mean_density + se_density), 
                position = position_dodge(width = 0.3), width = 0) +
  geom_point(size = 2, position = position_dodge(width = 0.3), alpha = 0.8) +
  labs(x = "Year", y = "Density") +
  facet_wrap(~region, scales = "free_y") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 
