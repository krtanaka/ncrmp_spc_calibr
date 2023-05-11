library(readr)
library(lubridate)
library(colorRamps)
library(ggplot2)
library(dplyr)

rm(list = ls())

region = c("MHI", "MARIAN", "NWHI", "PRIAs", "SAMOA")
var = c("abund", "biom")[1]

species = "APVI"
species = "ACLI"
species = "CAME"
species = "MOGR"
species = "NALI"
species = "SCSC"
species = "LUFU"

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

df_calibr <- df %>%
  group_by(island, depth, date_, latitude, longitude, species, region, year, month, day) %>%
  summarize(density = mean(density, na.rm = T)) %>% 
  mutate(method = "nSPC_BLT_TOW")

df = rbind(df, df_calibr)

save(df, file = paste0("output/calibr_df/calibr_", species, "_", var, ".RData"))
load(paste0("output/calibr_df/calibr_", species, "_", var, ".RData"))

if(var == "abund") unit = expression("Individuals (n) per 100" ~ m^2~"")
if(var == "biom") unit = expression("Biomass (g) per 100" ~ m^2~"")

df %>% 
  mutate(longitude = round(longitude, 1), 
         latitude = round(latitude, 1)) %>% 
  group_by(method, longitude, latitude) %>%
  summarise(density = mean(density)) %>%
  ggplot(aes(longitude, latitude)) + 
  geom_point(aes(size = density, fill = density), shape = 21, alpha = 0.5) +
  scale_fill_gradientn(colours = matlab.like(100), guide = "legend", trans = "sqrt") +
  facet_grid(~ method) +
  ggtitle(paste0(species, ": ", min(df$year), "-", max(df$year))) + 
  # labs(x = expression(paste("Longitude ", degree, "W", sep = "")),
  # y = expression(paste("Latitude ", degree, "N", sep = ""))) +
  guides(color = guide_legend(unit), 
         fill = guide_legend(unit),
         size = guide_legend(unit)) + 
  theme(legend.position = "bottom",
        legend.key = element_rect(colour = NA, fill = NA),
        legend.background = element_rect(fill = "transparent", colour = NA),
        legend.box.background = element_rect(fill = "transparent", colour = NA))

ggsave(last_plot(),file = paste0("output/plot/calibr_", species, "_map_a_", var, ".pdf"), height = 5, width = 12)

df %>% 
  filter(method == "nSPC_BLT_TOW") %>%
  mutate(longitude = round(longitude, 1), 
         latitude = round(latitude, 1)) %>% 
  group_by(longitude, latitude, region) %>%
  summarise(density = mean(density)) %>%
  ggplot(aes(longitude, latitude)) + 
  geom_point(aes(size = density, fill = density), shape = 21, alpha = 0.7) +
  scale_fill_gradientn(colours = matlab.like(100), guide = "legend", trans = "sqrt") +
  ggtitle(paste0(species, ": ", min(df$year), "-", max(df$year))) + 
  guides(color = guide_legend(unit), 
         fill = guide_legend(unit),
         size = guide_legend(unit)) + 
  coord_fixed() + 
  theme(legend.position = c(0.18, 0.2),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.background = element_rect(fill = "transparent", colour = NA),
        legend.box.background = element_rect(fill = "transparent", colour = NA))

ggsave(last_plot(),file = paste0("output/plot/calibr_", species, "_map_b_", var, ".pdf"), height = 5, width = 7)

df %>% 
  filter(region %in% c("MHI")) %>% 
  mutate(longitude = round(longitude, 1), 
         latitude = round(latitude, 1)) %>% 
  group_by(method, longitude, latitude, year) %>%
  summarise(density = mean(density)) %>%
  ggplot(aes(longitude, latitude)) + 
  geom_point(aes(size = density, fill = density), shape = 21, alpha = 0.8) +
  scale_fill_gradientn(colours = matlab.like(100), guide = "legend", trans = "sqrt") +
  facet_grid(method ~ year) +
  ggtitle(species) + 
  guides(color = guide_legend(unit), 
         fill = guide_legend(unit),
         size = guide_legend(unit)) + 
  theme(legend.position = "bottom",
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.background = element_rect(fill = "transparent", colour = NA),
        legend.box.background = element_rect(fill = "transparent", colour = NA))

ggsave(last_plot(),file = paste0("output/plot/calibr_", species, "_map_c_", var, ".pdf"),height = 7, width = 12)

df %>% 
  filter(method == "nSPC_BLT_TOW") %>%
  mutate(depth = round(depth, 1)) %>%
  group_by(method, region, depth) %>%
  summarise(density = mean(density, na.rm = T)) %>%
  ggplot(aes(depth, density)) + 
  geom_smooth(method = "gam", color = "gray60", fill = "gray80") +
  geom_point(aes(fill = density), shape = 21, alpha = 0.8, size = 3, show.legend = F) +
  scale_fill_gradientn(colours = matlab.like(100), trans = "sqrt") +
  labs(x = "Depth (m)", y = unit) +
  facet_grid(~region) +
  ggtitle(species) + 
  guides(color = guide_legend(unit), 
         fill = guide_legend(unit),
         size = guide_legend(unit))

ggsave(last_plot(),file = paste0("output/plot/calibr_", species, "_depth_", var, ".pdf"), height = 5, width = 10)

df %>%
  # filter(method != "nSPC_BLT_TOW") %>%
  mutate(YEAR = format(date_, "%Y")) %>% 
  group_by(year, method, region) %>%
  summarise(mean_density = mean(density), se_density = sd(density)/sqrt(n())) %>%
  mutate(mean_density = ifelse(mean_density == 0, NA, mean_density),
         se_density = ifelse(se_density == 0, NA, se_density)) %>% 
  ggplot(aes(x = year, y = mean_density, color = method)) +
  geom_errorbar(aes(ymin = mean_density - se_density, ymax = mean_density + se_density), 
                position = position_dodge(width = 0.5), width = 0, show.legend = F) +
  geom_point(size = 2, position = position_dodge(width = 0.5)) +
  scale_color_discrete("") + 
  ggtitle(species) + 
  labs(x = NULL, y = unit) +
  facet_wrap(~region, scales = "free_y") +
  scale_x_discrete(limits = unique(df$year)) + # Add this line
  theme(legend.position = c(0.85, 0.25), 
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.key = element_rect(colour = NA, fill = NA),
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.box.background = element_rect(fill = "transparent", colour = NA))

ggsave(last_plot(),file = paste0("output/plot/calibr_", species, "_ts_a_", var, ".pdf"), height = 5, width = 10)

df %>%
  filter(method != "nSPC_BLT_TOW") %>%
  mutate(YEAR = format(date_, "%Y")) %>% 
  group_by(year, region) %>%
  summarise(mean_density = mean(density), se_density = sd(density)/sqrt(n())) %>%
  mutate(mean_density = ifelse(mean_density == 0, NA, mean_density),
         se_density = ifelse(se_density == 0, NA, se_density)) %>% 
  ggplot(aes(x = year, y = mean_density, fill = mean_density)) +
  geom_errorbar(aes(ymin = mean_density - se_density, ymax = mean_density + se_density), width = 0, show.legend = F) +
  geom_point(size = 3, shape = 21, show.legend = F) +
  scale_fill_gradientn(colours = matlab.like(100), "", tran = "sqrt") +
  ggtitle(species) + 
  labs(x = NULL, y = unit) +
  facet_wrap(~region, scales = "free_y") +
  scale_x_discrete(limits = unique(df$year)) + # Add this line
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave(last_plot(),file = paste0("output/plot/calibr_", species, "_ts_b_", var, ".pdf"), height = 5, width = 10)

df %>%
  filter(method != "nSPC_BLT_TOW") %>% 
  mutate(YEAR = format(date_, "%Y")) %>% 
  group_by(year, island) %>%
  summarise(mean_density = mean(density), se_density = sd(density)/sqrt(n())) %>%
  mutate(mean_density = ifelse(mean_density == 0, NA, mean_density),
         se_density = ifelse(se_density == 0, NA, se_density)) %>% 
  na.omit() %>% 
  ggplot(aes(x = year, y = mean_density, fill = mean_density)) +
  geom_errorbar(aes(ymin = mean_density - se_density, ymax = mean_density + se_density), 
                width = 0, position = position_dodge(width = 0.5), show.legend = F) +
  geom_point(size = 3, shape = 21, position = position_dodge(width = 0.5), show.legend = F) +
  scale_fill_gradientn(colours = matlab.like(100), "", tran = "sqrt") +
  ggtitle(species) + 
  labs(x = NULL, y = unit) + 
  facet_wrap(~island, scales = "free_y") +
  scale_x_discrete(limits = unique(df$year)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave(last_plot(),file = paste0("output/plot/calibr_", species, "_ts_c_", var, ".pdf"), height = 10, width = 20)
