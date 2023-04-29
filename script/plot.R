library(readr)
library(lubridate)
library(colorRamps)
library(ggplot2)

rm(list = ls())

region = c("MHI", "MARIAN", "NWHI", "PRIAs", "SAMOA")[1]
var = c("abund", "biom")[1]
species = "APVI"

calibr_belt = read_csv(paste0("output/spc_belt_", var, "_", region, "_GLMM/summary_table.csv")) %>% 
  subset(GROUP == species) %>% 
  subset(METHOD != "1_nSPC")

calibr_tow = read_csv(paste0("output/spc_tow_", var, "_", region, "_GLMM/summary_table.csv")) %>% 
  subset(GROUP == species) %>% 
  subset(METHOD != "1_nSPC")

calibr_tow = calibr_tow$POS.GCF
calibr_belt = calibr_belt$POS.GCF

belt <- readRDS(paste0("data/belt.site.", var, ".size.20002009.", region, ".rds")) %>% 
  group_by(ISLAND, DEPTH, METHOD, DATE_, LATITUDE, LONGITUDE, SPECIES) %>% 
  summarise(DENSITY = sum(!!sym(paste0(var, ".site"))), .groups = "drop")  %>% 
  subset(SPECIES == species) %>%
  mutate(DENSITY = (DENSITY / calibr_belt))

spc = readRDS(paste0("data/nSPC.site.", var, ".size.20092022.", region, ".rds"))  %>% 
  group_by(ISLAND, DEPTH, METHOD, DATE_, LATITUDE, LONGITUDE, SPECIES) %>% 
  summarise(DENSITY = sum(!!sym(paste0(var, ".site"))), .groups = "drop") %>% 
  subset(SPECIES == species)

tow = readRDS(paste0("data/tow.segment.", var, ".size.20002017.", region, ".rds"))  %>% 
  subset(CENTROIDLON != 0) %>%
  mutate(LONGITUDE = CENTROIDLON,
         LATITUDE = CENTROIDLAT) %>% 
  group_by(ISLAND, DEPTH, METHOD, DATE_, LATITUDE, LONGITUDE, SPECIES) %>% 
  summarise(DENSITY = sum(!!sym(paste0(var, ".segment"))), .groups = "drop") %>% 
  subset(SPECIES == species) %>%
  mutate(DENSITY = (DENSITY / calibr_tow))

df = rbind(spc, belt, tow)

df$YEAR <- year(df$DATE_)
df$MONTH <- month(df$DATE_)

df %>% 
  mutate(DENSITY = DENSITY*100,
         LONGITUDE = round(LONGITUDE, 1), 
         LATITUDE= round(LATITUDE, 1)) %>% 
  group_by(METHOD, LONGITUDE, LATITUDE) %>%
  summarise(DENSITY = mean(DENSITY)) %>%
  ggplot(aes(LONGITUDE, LATITUDE)) + 
  geom_point(aes(size = DENSITY, fill = DENSITY, color = DENSITY), shape = 21, alpha = 0.7) +
  scale_color_gradientn(colours = matlab.like(100), guide = "legend", trans = "sqrt") +
  scale_fill_gradientn(colours = matlab.like(100), guide = "legend", trans = "sqrt") +
  # facet_grid(METHOD~YEAR) +
  facet_grid(~METHOD) +
  labs(x = expression(paste("Longitude ", degree, "W", sep = "")),
       y = expression(paste("Latitude ", degree, "N", sep = ""))) +
  guides(color = guide_legend(expression("Individuals per 100" ~ m^2~"")), 
         fill = guide_legend(expression("Individuals per 100" ~ m^2~"")),
         size = guide_legend(expression("Individuals per 100" ~ m^2~""))) + 
  theme(legend.position = "bottom",
        panel.background = element_rect(fill = "gray10"),
        panel.grid = element_line(color = "gray20", linewidth = 0.1)) 

df %>% 
  mutate(DENSITY = DENSITY*100,
         LONGITUDE = round(LONGITUDE, 1), 
         LATITUDE= round(LATITUDE, 1)) %>% 
  group_by(METHOD, LONGITUDE, LATITUDE, YEAR) %>%
  summarise(DENSITY = mean(DENSITY)) %>%
  ggplot(aes(LONGITUDE, LATITUDE)) + 
  geom_point(aes(size = DENSITY, fill = DENSITY, color = DENSITY), shape = 21, alpha = 0.7) +
  scale_color_gradientn(colours = matlab.like(100), guide = "legend", trans = "sqrt") +
  scale_fill_gradientn(colours = matlab.like(100), guide = "legend", trans = "sqrt") +
  facet_grid(METHOD~YEAR) +
  labs(x = expression(paste("Longitude ", degree, "W", sep = "")),
       y = expression(paste("Latitude ", degree, "N", sep = ""))) +
  guides(color = guide_legend(expression("Individuals per 100" ~ m^2~"")), 
         fill = guide_legend(expression("Individuals per 100" ~ m^2~"")),
         size = guide_legend(expression("Individuals per 100" ~ m^2~""))) + 
  theme(legend.position = "bottom",
        panel.background = element_rect(fill = "gray10"),
        panel.grid = element_line(color = "gray20", linewidth = 0.1)) 

df %>% 
  mutate(DENSITY = DENSITY*100) %>% 
  ggplot(aes(DEPTH, DENSITY)) + 
  geom_point(aes(size = DENSITY, fill = DENSITY, color = DENSITY), shape = 21, alpha = 0.7) +
  scale_color_gradientn(colours = matlab.like(100), guide = "legend", trans = "sqrt") +
  scale_fill_gradientn(colours = matlab.like(100), guide = "legend", trans = "sqrt") +
  facet_grid(~METHOD) +
  guides(color = guide_legend(expression("Individuals per 100" ~ m^2~"")), 
         fill = guide_legend(expression("Individuals per 100" ~ m^2~"")),
         size = guide_legend(expression("Individuals per 100" ~ m^2~""))) + 
  theme(legend.position = "bottom",
        panel.background = element_rect(fill = "gray10"),
        panel.grid = element_line(color = "gray20")) 

df %>%
  mutate(YEAR = format(DATE_, "%Y")) %>%
  group_by(YEAR, METHOD) %>%
  summarise(mean_density = mean(DENSITY), se_density = sd(DENSITY)/sqrt(n())) %>%
  mutate(mean_density = ifelse(mean_density == 0, NA, mean_density),
         se_density = ifelse(se_density == 0, NA, se_density)) %>% 
  ggplot(aes(x = YEAR, y = mean_density, color = METHOD, group = METHOD)) +
  geom_errorbar(aes(ymin = mean_density - se_density, ymax = mean_density + se_density), position=position_dodge(width = 0.5), width = 0) +
  geom_point(size = 3, position=position_dodge(width=0.5)) +
  labs(x = "Year", y = "Density") +
  theme(legend.position = "bottom",
        panel.background = element_rect(fill = "gray10"),
        panel.grid = element_line(color = "gray20")) 
