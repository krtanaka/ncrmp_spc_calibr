library(readr)
library(lubridate)
library(colorRamps)
library(ggplot2)
library(dplyr)
library(boot)
library(ggh4x)

rm(list = ls())

var = c("abund", "biom")[1]

region = c("MHI", "MARIAN", "NWHI", "PRIAs", "SAMOA")[1]

species = c(
  "APXA",
  "APVI",
  "ACLI",
  "ACNC",
  "CAME",
  "MOGR",
  "NALI",
  "SCSC",
  "LUFU",
  "LUKA",
  "CHUN",
  "BOMU")

for (s in 1:length(species)) {
  
  # s = 1
  
  df = NULL
  
  for (r in 1:length(region)) {
    
    # r = 1
    
    calibr_belt = read_csv(paste0("output/spc_belt_", var, "_", region[r], "_GLMM/summary_table.csv")) %>% 
      filter(GROUP == species[s]) %>%
      filter(METHOD != "1_nSPC")
    
    calibr_tow = read_csv(paste0("output/spc_tow_", var, "_", region[r], "_GLMM/summary_table.csv")) %>% 
      filter(GROUP == species[s]) %>%
      filter(METHOD != "1_nSPC")
    
    calibr_tow_pres = calibr_tow$PRES.GCF
    calibr_belt_pres = calibr_belt$PRES.GCF
    
    calibr_tow_pos = calibr_tow$POS.GCF
    calibr_belt_pos = calibr_belt$POS.GCF
    
    belt <- readRDS(paste0("data/belt.site.", var, ".size.20002009.", region[r], ".rds")) %>%
      filter(SPECIES == species[s]) %>%
      group_by(ISLAND, DEPTH, METHOD, DATE_, LATITUDE, LONGITUDE, SPECIES) %>%
      summarise(DENSITY = sum(!!sym(paste0(var, ".site"))), .groups = "drop") %>%
      mutate(PRESENCE = inv.logit(logit(DENSITY > 0) - calibr_belt_pres),
             DENSITY = DENSITY / calibr_belt_pos)
    
    if (dim(belt)[1] == 0) {
      
      belt <- readRDS(paste0("data/belt.site.", var, ".size.20002009.", region[r], ".rds")) %>%
        group_by(ISLAND, DEPTH, METHOD, DATE_, LATITUDE, LONGITUDE) %>%
        summarise(SPECIES = species[s]) %>% 
        mutate(DENSITY = 0,
               PRESENCE = 0)
      
    }
    
    tow = readRDS(paste0("data/tow.segment.", var, ".size.20002017.", region[r], ".rds"))  %>% 
      filter(SPECIES == species[s] & CENTROIDLON != 0 & SIZE_10cm != "(40,50]") %>%
      mutate(LONGITUDE = CENTROIDLON,
             LATITUDE = CENTROIDLAT) %>% 
      group_by(ISLAND, DEPTH, METHOD, DATE_, LATITUDE, LONGITUDE, SPECIES) %>% 
      summarise(DENSITY = sum(!!sym(paste0(var, ".segment"))), .groups = "drop") %>% 
      mutate(PRESENCE = inv.logit(logit(DENSITY > 0) - calibr_tow_pres),
             DENSITY = DENSITY / calibr_tow_pos)
    
    if (dim(tow)[1] == 0) {
      
      tow = readRDS(paste0("data/tow.segment.", var, ".size.20002017.", region[r], ".rds"))  %>% 
        filter(CENTROIDLON != 0 & SIZE_10cm != "(40,50]") %>%
        mutate(LONGITUDE = CENTROIDLON,
               LATITUDE = CENTROIDLAT) %>% 
        group_by(ISLAND, DEPTH, METHOD, DATE_, LATITUDE, LONGITUDE) %>% 
        summarise(SPECIES = species[s]) %>% 
        mutate(DENSITY = 0,
               PRESENCE = 0)
      
    }
    
    spc = readRDS(paste0("data/nSPC.site.", var, ".size.20092022.", region[r], ".rds"))  %>% 
      filter(SPECIES == species[s]) %>%
      group_by(ISLAND, DEPTH, METHOD, DATE_, LATITUDE, LONGITUDE, SPECIES) %>% 
      summarise(DENSITY = sum(!!sym(paste0(var, ".site"))), .groups = "drop") %>%
      mutate(PRESENCE = as.integer(DENSITY > 0)) 
    
    if (dim(spc)[1] == 0) {
      
      spc = readRDS(paste0("data/nSPC.site.", var, ".size.20092022.", region[r], ".rds"))  %>% 
        group_by(ISLAND, DEPTH, METHOD, DATE_, LATITUDE, LONGITUDE) %>% 
        summarise(SPECIES = species[s]) %>% 
        mutate(DENSITY = 0,
               PRESENCE = 0)
      
    }
    
    df_i = rbind(spc, belt, tow)
    df_i$region = region[r]
    
    df = rbind(df, df_i)
    
  }
  
  df = df %>% rename_all(tolower)
  
  df_tow_blt <- df %>%
    filter(method != "nSPC") %>% 
    mutate(density = density * presence)
  
  df_spc = df %>%
    filter(method == "nSPC")
  
  df_nSPC_BLT_TOW = rbind(df_spc, df_tow_blt) %>% 
    group_by(region, island, depth, date_, latitude, longitude, species) %>%
    summarize(presence = mean(presence, na.rm = T),
              density = mean(density, na.rm = T)) %>% 
    mutate(method = "nSPC_BLT_TOW")
  
  df <- rbind(df_spc, df_tow_blt, df_nSPC_BLT_TOW) %>% 
    mutate(year = year(date_),
           month = month(date_),
           day = day(date_),
           longitude = ifelse(longitude < 0, longitude + 360, longitude),
           density = density * 100)
  
  save(df, file = paste0("output/calibr_df/calibr_", species[s], "_", var, ".RData"))
  
  if(var == "abund") unit = expression("Individuals (n) per 2.4" ~ km^2~"")
  if(var == "biom") unit = expression("Biomass (g) per 2.4" ~ km^2~"")
  
  df %>% 
    # filter(density > 0) %>%
    filter(method != "nSPC_BLT_TOW") %>%
    mutate(longitude = round(longitude, 1),
           latitude = round(latitude, 1)) %>%
    # mutate(longitude = round(longitude / 0.5) * 0.5,
           # latitude = round(latitude / 0.5) * 0.5) %>%
    group_by(method, island, longitude, latitude) %>%
    summarise(density = mean(density)) %>%
    ggplot(aes(longitude, latitude)) + 
    geom_polygon(data = fortify(maps::map("world2", plot = F, fill = T)),
                 aes(x = ifelse(long < 0, long + 360, long), y = lat, group = group)) +
    geom_point(aes(size = density, fill = density), shape = 21, alpha = 0.5) +
    scale_fill_gradientn(colours = matlab.like(100), guide = "legend") +
    facet_wrap(~ method) +
    ggtitle(paste0(species[s], ": ", min(df$year), "-", max(df$year)),
            subtitle = paste0("Maximum observation per 100 sq.m = ", round(max(df$density), 1))) + 
    coord_equal(xlim = range(df$longitude), ylim = range(df$latitude)) +
    labs(x = expression(paste("Longitude ", degree, "W", sep = "")),
         y = expression(paste("Latitude ", degree, "N", sep = ""))) +
    guides(color = guide_legend(unit), 
           fill = guide_legend(unit),
           size = guide_legend(unit)) + 
    theme(legend.position = "bottom",
          legend.key = element_rect(colour = NA, fill = NA),
          legend.background = element_rect(fill = "transparent", colour = NA),
          legend.box.background = element_rect(fill = "transparent", colour = NA))
  
  # ggsave(last_plot(),file = paste0("output/plot/map_a_", species[s], "_", var, ".pdf"), height = 8, width = 16)
  ggsave(last_plot(),file = paste0("output/plot/map_a_", species[s], "_", var, ".png"), height = 6, width = 18, units = "in")
  
  df %>% 
    filter(density > 0) %>%
    filter(method == "nSPC_BLT_TOW") %>%
    mutate(longitude = round(longitude, 1),
           latitude = round(latitude, 1)) %>%
    # mutate(longitude = round(longitude / 0.5) * 0.5,
    # latitude = round(latitude / 0.5) * 0.5) %>%
    group_by(longitude, latitude, region) %>%
    summarise(density = mean(density)) %>%
    ggplot(aes(longitude, latitude)) + 
    geom_polygon(data = fortify(maps::map("world2", plot = F, fill = T)), aes(x = ifelse(long < 0, long + 360, long), y = lat, group = group)) +
    geom_point(aes(size = density, fill = density), shape = 21, alpha = 0.5) +
    scale_fill_gradientn(colours = matlab.like(100), guide = "legend", trans = "sqrt") +
    ggtitle(paste0(species[s], ": ", min(df$year), "-", max(df$year)),
            subtitle = paste0("Maximum observation per 100 sq.m = ", round(max(df$density), 1))) + 
    coord_equal(xlim = range(df$longitude), ylim = range(df$latitude)) +
    labs(x = expression(paste("Longitude ", degree, "", sep = "")),
         y = expression(paste("Latitude ", degree, "", sep = ""))) +
    guides(color = guide_legend(unit), 
           fill = guide_legend(unit),
           size = guide_legend(unit)) + 
    theme(legend.position = c(0.15, 0.2),
          legend.key = element_rect(colour = NA, fill = NA),
          legend.background = element_rect(fill = "transparent", colour = NA),
          legend.box.background = element_rect(fill = "transparent", colour = NA))
  
  # ggsave(last_plot(),file = paste0("output/plot/map_b_", species[s], "_", var, ".pdf"), height = 8, width = 16)
  ggsave(last_plot(),file = paste0("output/plot/map_b_", species[s], "_", var, ".png"), height = 6, width = 8, units = "in")
  
  df %>% 
    # filter(density > 0) %>%
    mutate(longitude = round(longitude, 1),
           latitude = round(latitude, 1)) %>%
    # mutate(longitude = round(longitude / 0.5) * 0.5,
    # latitude = round(latitude / 0.5) * 0.5) %>%
    group_by(method, longitude, latitude, year) %>%
    summarise(density = mean(density)) %>%
    ggplot(aes(longitude, latitude)) + 
    geom_point(aes(size = density, fill = density), shape = 21, alpha = 0.5) +
    scale_fill_gradientn(colours = matlab.like(100), guide = "legend", trans = "sqrt") +
    facet_grid(method ~ year) +
    ggtitle(paste0(species[s], ": ", min(df$year), "-", max(df$year)),
            subtitle = paste0("Maximum density per 100 sq.m : ", round(max(df$density), 1))) + 
    guides(color = guide_legend(unit), 
           fill = guide_legend(unit),
           size = guide_legend(unit)) + 
    theme(legend.position = "bottom",
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          legend.key = element_rect(colour = NA, fill = NA),
          legend.background = element_rect(fill = "transparent", colour = NA),
          legend.box.background = element_rect(fill = "transparent", colour = NA))
  
  # ggsave(last_plot(),file = paste0("output/plot/map_c_", species[s], "_", var, ".pdf"),height = 10, width = 30)
  ggsave(last_plot(),file = paste0("output/plot/map_c_", species[s], "_", var, ".png"),height = 7, width = 28, units = "in")
  
  if(var == "abund") unit = expression("Individuals (n) per 100" ~ m^2~"")
  if(var == "biom") unit = expression("Biomass (g) per 100" ~ m^2~"")
  
  df %>%
    # filter(density > 0) %>%
    filter(method == "nSPC_BLT_TOW") %>%
    # mutate(depth = round(depth, 1)) %>%
    # group_by(method, region, depth) %>%
    # summarise(density = mean(density, na.rm = T)) %>%
    ggplot(aes(depth, density)) + 
    # geom_smooth(method = "gam", color = "gray60", fill = "gray80") +
    geom_point(aes(fill = density), shape = 21, alpha = 0.8, size = 3, show.legend = F) +
    scale_fill_gradientn(colours = matlab.like(100), trans = "sqrt") +
    labs(x = "Depth (m)", y = unit) +
    facet_grid(~region) +
    ggtitle(species[s]) + 
    guides(color = guide_legend(unit), 
           fill = guide_legend(unit),
           size = guide_legend(unit))
  
  # ggsave(last_plot(),file = paste0("output/plot/depth_", species[s], "_", var, ".pdf"), height = 5, width = 10)
  ggsave(last_plot(),file = paste0("output/plot/depth_", species[s], "_", var, ".png"), height = 5, width = 10, units = "in")
  
  df %>%
    filter(method != "nSPC_BLT_TOW") %>%
    mutate(YEAR = format(date_, "%Y")) %>% 
    group_by(year, method, region) %>%
    summarise(mean_density = mean(density), se_density = sd(density)/sqrt(n())) %>%
    # mutate(mean_density = ifelse(mean_density == 0, NA, mean_density),
    #        se_density = ifelse(se_density == 0, NA, se_density)) %>% 
    ggplot(aes(x = year, y = mean_density, color = method, group = method)) +
    geom_errorbar(aes(ymin = mean_density - se_density, ymax = mean_density + se_density), 
                  position = position_dodge(width = 0.5), width = 0, show.legend = F) +
    geom_point(size = 2, position = position_dodge(width = 0.5)) +
    geom_line(position = position_dodge(width = 0.5)) + 
    scale_color_discrete("") + 
    ggtitle(species[s]) + 
    labs(x = NULL, y = unit) +
    facet_wrap(~region, scales = "free_y") +
    scale_x_discrete(limits = unique(df$year)) + # Add this line
    theme(legend.position = c(0.85, 0.25), 
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.key = element_rect(colour = NA, fill = NA),
          legend.background = element_rect(fill = "transparent", colour = NA),
          legend.box.background = element_rect(fill = "transparent", colour = NA))
  
  # ggsave(last_plot(),file = paste0("output/plot/ts_a_", species[s], "_", var, ".pdf"), height = 5, width = 10)
  ggsave(last_plot(),file = paste0("output/plot/ts_a_", species[s], "_", var, ".png"), height = 5, width = 10, units = "in")
  
  df %>%
    filter(method != "nSPC_BLT_TOW") %>%
    mutate(YEAR = format(date_, "%Y")) %>% 
    group_by(year, region) %>%
    summarise(mean_density = mean(density), se_density = sd(density)/sqrt(n())) %>%
    # mutate(mean_density = ifelse(mean_density == 0, NA, mean_density),
    # se_density = ifelse(se_density == 0, NA, se_density)) %>%
    ggplot(aes(x = year, y = mean_density, fill = mean_density)) +
    geom_errorbar(aes(ymin = mean_density - se_density, ymax = mean_density + se_density), width = 0, show.legend = F) +
    geom_point(size = 3, shape = 21, show.legend = F) +
    scale_fill_gradientn(colours = matlab.like(100), "", tran = "sqrt") +
    ggtitle(species[s]) + 
    labs(x = NULL, y = unit) +
    facet_wrap(~region, scales = "free_y") +
    scale_x_discrete(limits = unique(df$year)) + # Add this line
    theme(legend.position = "bottom", 
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  # ggsave(last_plot(),file = paste0("output/plot/ts_b_", species[s], "_", var, ".pdf"), height = 5, width = 10)
  ggsave(last_plot(),file = paste0("output/plot/ts_b_", species[s], "_", var, ".png"), height = 5, width = 10, units = "in")
  
  dfi = df %>%
    filter(method != "nSPC_BLT_TOW") %>% 
    mutate(YEAR = format(date_, "%Y")) %>% 
    group_by(year, island) %>%
    summarise(mean_density = mean(density), se_density = sd(density)/sqrt(n())) %>%
    # mutate(mean_density = ifelse(mean_density == 0, NA, mean_density),
    # se_density = ifelse(se_density == 0, NA, se_density)) %>%
    na.omit() 
  
  dfi %>%
    ggplot(aes(x = year, y = mean_density, fill = mean_density)) +
    geom_errorbar(aes(ymin = mean_density - se_density, ymax = mean_density + se_density), 
                  width = 0, position = position_dodge(width = 0.5), show.legend = F) +
    geom_point(size = 3, shape = 21, position = position_dodge(width = 0.5), show.legend = F) +
    scale_fill_gradientn(colours = matlab.like(100), "", tran = "sqrt") +
    ggtitle(species[s]) + 
    labs(x = NULL, y = unit) + 
    facet_wrap(~island, scales = "free_y") +
    scale_x_discrete(limits = unique(dfi$year)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  # ggsave(last_plot(),file = paste0("output/plot/ts_c_", species[s], "_", var, ".pdf"), height = 5, width = 10)
  ggsave(last_plot(),file = paste0("output/plot/ts_c_", species[s], "_", var, ".png"), height = 10, width = 20, units = "in")
  
}
