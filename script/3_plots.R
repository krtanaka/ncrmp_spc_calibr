library(readr)
library(lubridate)
library(colorRamps)
library(ggplot2)
library(dplyr)
library(boot)
library(ggh4x)
library(classInt)
library(viridis)
library(scales)

rm(list = ls())

var = c("abund", "biom")[1]

region = c("MHI", "MARIAN", "NWHI", "PRIAs", "SAMOA")[1]

species = c(
  "ACLI", #1
  "ACNC", #2
  "APVI", #3
  "APXA", #4
  "CAIG", #5
  "CALU", #6
  "CAME", #7
  "CEAR", #8
  "LERU", #9
  "LUBO", #10
  "LUFU", #11
  "LUKA", #12
  "MOGR", #13
  "NALI", #14
  "NAUN", #15
  "SCLY", #16
  "SCSC", #17
  "VALO"  #18
)[3]

# -------------------------------------------------------------------------
# Plot helpers
# -------------------------------------------------------------------------

region_labs <- c(
  MHI = "Main Hawaiian Islands",
  MARIAN = "Mariana Archipelago",
  NWHI = "Northwestern Hawaiian Islands",
  PRIAs = "Pacific Remote Island Areas",
  SAMOA = "American Samoa"
)

taxa_labs <- c(
  APVI = "Aprion virescens (Uku)"
)

method_labs <- c(
  nSPC = "Stationary Point Count",
  BLT = "Belt transect",
  TOW = "Towed diver",
  nSPC_BLT_TOW = "Calibrated composite"
)

theme_pub <- theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", size = 18, margin = margin(b = 4)),
    plot.subtitle = element_text(size = 11.5, colour = "grey30", lineheight = 1.15, margin = margin(b = 10)),
    plot.caption = element_text(size = 9.5, colour = "grey40", hjust = 0, margin = margin(t = 10)),
    axis.title = element_text(size = 12.5),
    axis.text = element_text(size = 10.5, colour = "grey30"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey86", linewidth = 0.35),
    strip.background = element_rect(fill = "grey94", colour = NA),
    strip.text = element_text(face = "bold"),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    legend.key = element_blank(),
    legend.box.margin = margin(0, 0, 0, 8),
    plot.margin = margin(12, 16, 12, 12)
  )

prep_map_bins <- function(pdat) {
  cap <- quantile(pdat$density, 0.99, na.rm = TRUE)
  pdat <- pdat %>% mutate(density_c = pmin(density, cap))
  brks <- signif(quantile(pdat$density_c, c(.10, .25, .50, .75, .90, .98), na.rm = TRUE), 2)
  brks <- sort(unique(c(0, brks, cap)))
  if (length(brks) < 3) brks <- sort(unique(c(0, pretty(c(0, cap), n = 4), cap)))
  
  bin_labs <- paste0(
    scales::number(brks[-length(brks)], accuracy = 0.01), "–",
    scales::number(brks[-1], accuracy = 0.01)
  )
  
  pdat$density_bin <- cut(pdat$density_c, breaks = brks, include.lowest = TRUE, labels = bin_labs)
  lv <- levels(pdat$density_bin)
  
  list(
    dat = pdat,
    levels = lv,
    fill_vals = setNames(viridis::magma(length(lv), direction = -1), lv),
    size_vals = setNames(seq(2.2, 9.2, length.out = length(lv)), lv),
    cap = cap
  )
}

facet_dims <- function(n_panels, ncol, panel_w = 3.7, panel_h = 2.8) {
  nrow <- ceiling(n_panels / ncol)
  c(width = max(6, 2 + panel_w * ncol), height = max(4.8, 1.8 + panel_h * nrow))
}

for (s in 1:length(species)) {
  
  df = NULL
  
  for (r in 1:length(region)) {
    
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
        mutate(DENSITY = 0, PRESENCE = 0)
    }
    
    tow = readRDS(paste0("data/tow.segment.", var, ".size.20002017.", region[r], ".rds")) %>%
      filter(SPECIES == species[s] & CENTROIDLON != 0 & SIZE_10cm != "(40,50]") %>%
      mutate(LONGITUDE = CENTROIDLON, LATITUDE = CENTROIDLAT) %>%
      group_by(ISLAND, DEPTH, METHOD, DATE_, LATITUDE, LONGITUDE, SPECIES) %>%
      summarise(DENSITY = sum(!!sym(paste0(var, ".segment"))), .groups = "drop") %>%
      mutate(PRESENCE = inv.logit(logit(DENSITY > 0) - calibr_tow_pres),
             DENSITY = DENSITY / calibr_tow_pos)
    
    if (dim(tow)[1] == 0) {
      tow = readRDS(paste0("data/tow.segment.", var, ".size.20002017.", region[r], ".rds")) %>%
        filter(CENTROIDLON != 0 & SIZE_10cm != "(40,50]") %>%
        mutate(LONGITUDE = CENTROIDLON, LATITUDE = CENTROIDLAT) %>%
        group_by(ISLAND, DEPTH, METHOD, DATE_, LATITUDE, LONGITUDE) %>%
        summarise(SPECIES = species[s]) %>%
        mutate(DENSITY = 0, PRESENCE = 0)
    }
    
    spc = readRDS(paste0("data/nSPC.site.", var, ".size.20092022.", region[r], ".rds")) %>%
      filter(SPECIES == species[s]) %>%
      group_by(ISLAND, DEPTH, METHOD, DATE_, LATITUDE, LONGITUDE, SPECIES) %>%
      summarise(DENSITY = sum(!!sym(paste0(var, ".site"))), .groups = "drop") %>%
      mutate(PRESENCE = as.integer(DENSITY > 0))
    
    if (dim(spc)[1] == 0) {
      spc = readRDS(paste0("data/nSPC.site.", var, ".size.20092022.", region[r], ".rds")) %>%
        group_by(ISLAND, DEPTH, METHOD, DATE_, LATITUDE, LONGITUDE) %>%
        summarise(SPECIES = species[s]) %>%
        mutate(DENSITY = 0, PRESENCE = 0)
    }
    
    df_i = rbind(spc, belt, tow)
    df_i$region = region[r]
    df = rbind(df, df_i)
  }
  
  df = df %>% rename_all(tolower)
  
  df_tow_blt <- df %>%
    filter(method != "nSPC") %>%
    mutate(density = density * presence)
  
  df_spc = df %>% filter(method == "nSPC")
  
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
  
  save(df, file = paste0("output/calibr_df/calibr_", species[s], "_", var, "_", region, ".RData"))
  
  if (var == "abund") unit <- expression("Individuals (n) per 100" ~ m^2)
  if (var == "biom") unit <- expression("Biomass (g) per 100" ~ m^2)
  
  title_text <- if (species[s] %in% names(taxa_labs)) taxa_labs[[species[s]]] else species[s]
  region_text <- paste(unname(region_labs[unique(df$region)]), collapse = ", ")
  year_text <- paste0(min(df$year, na.rm = TRUE), "–", max(df$year, na.rm = TRUE))
  region_tag <- paste(unique(df$region), collapse = "-")
  
  max_obs_text <- function(x) {
    if (var == "abund")
      paste0("Maximum observed abundance = ", round(max(x, na.rm = TRUE), 1), " individuals per 100 m²")
    else
      paste0("Maximum observed biomass = ", round(max(x, na.rm = TRUE), 1), " g per 100 m²")
  }
  
  # -----------------------------------------------------------------------
  # MAP A: survey-specific spatial distributions
  # -----------------------------------------------------------------------
  
  pdat <- df %>%
    filter(density > 0, method != "nSPC_BLT_TOW") %>%
    mutate(longitude = round(longitude, 1), latitude = round(latitude, 1)) %>%
    group_by(method, island, longitude, latitude) %>%
    summarise(density = mean(density), .groups = "drop")
  
  pm <- prep_map_bins(pdat)
  pdat <- pm$dat
  lon_range <- range(pdat$longitude, na.rm = TRUE)
  lat_range <- range(pdat$latitude, na.rm = TRUE)
  lon_buffer <- diff(lon_range) * 0.06
  lat_buffer <- diff(lat_range) * 0.08
  n_methods <- length(unique(pdat$method))
  
  p = ggplot(pdat, aes(longitude, latitude)) +
    geom_polygon(
      data = fortify(maps::map("world2", plot = FALSE, fill = TRUE)),
      aes(x = ifelse(long < 0, long + 360, long), y = lat, group = group),
      fill = "grey25", colour = "white", linewidth = 0.25, inherit.aes = FALSE
    ) +
    geom_point(aes(size = density_bin, fill = density_bin),
               shape = 21, alpha = 0.85, colour = "grey20", stroke = 0.25) +
    scale_fill_manual(values = pm$fill_vals, breaks = pm$levels, name = unit, drop = FALSE) +
    scale_size_manual(values = pm$size_vals, breaks = pm$levels, name = unit, drop = FALSE) +
    guides(fill = guide_legend(title.position = "top",
                               override.aes = list(shape = 21, size = unname(pm$size_vals),
                                                   colour = "grey20", alpha = 0.9)),
           size = "none") +
    facet_wrap(~method, nrow = 1, labeller = labeller(method = method_labs)) +
    coord_equal(xlim = c(lon_range[1] - lon_buffer, lon_range[2] + lon_buffer),
                ylim = c(lat_range[1] - lat_buffer, lat_range[2] + lat_buffer),
                expand = FALSE) +
    labs(
      title = title_text,
      subtitle = paste0(region_text, " | Survey-specific calibrated ", ifelse(var == "abund", "abundance", "biomass"),
                        ", ", year_text, "\n", max_obs_text(pdat$density)),
      x = expression(Longitude ~ degree * W),
      y = expression(Latitude ~ degree * N),
      caption = "BLT and TOW observations standardized to the nSPC reference method. Mean values within 0.1° cells; values above the 99th percentile capped for display."
    ) +
    theme_pub +
    theme(legend.position = "right")
  
  dims <- facet_dims(n_methods, n_methods, 4.3, 4.1)
  
  ggsave(
    filename = paste0("output/plot/map_a_", species[s], "_", var, "_", region_tag, ".png"),
    plot = p, height = dims["height"], width = dims["width"],
    units = "in", dpi = 300, bg = "white"
  )
  
  pdat <- df %>%
    filter(density > 0, method == "nSPC_BLT_TOW") %>%
    mutate(longitude = round(longitude, 1), latitude = round(latitude, 1)) %>%
    group_by(island, longitude, latitude) %>%
    summarise(density = mean(density, na.rm = TRUE), .groups = "drop")
  
  pm <- prep_map_bins(pdat)
  pdat <- pm$dat
  lon_range <- range(pdat$longitude, na.rm = TRUE)
  lat_range <- range(pdat$latitude, na.rm = TRUE)
  lon_buffer <- diff(lon_range) * 0.06
  lat_buffer <- diff(lat_range) * 0.08
  
  p <- ggplot() +
    geom_polygon(
      data = fortify(maps::map("world2", plot = FALSE, fill = TRUE)),
      aes(x = ifelse(long < 0, long + 360, long), y = lat, group = group),
      fill = "grey25", colour = "white", linewidth = 0.25
    ) +
    geom_point(
      data = pdat,
      aes(x = longitude, y = latitude, size = density_bin, fill = density_bin),
      shape = 21, alpha = 0.8, colour = "grey20", stroke = 0.25
    ) +
    scale_fill_manual(values = pm$fill_vals, breaks = pm$levels, name = unit, drop = FALSE) +
    scale_size_manual(values = pm$size_vals, breaks = pm$levels, name = unit, drop = FALSE) +
    guides(
      fill = guide_legend(
        title.position = "top",
        override.aes = list(
          shape = 21, 
          size = unname(pm$size_vals),
          colour = "grey20", 
          alpha = 0.9
        )
      ),
      size = "none"
    ) +
    coord_equal(
      xlim = c(lon_range[1] - lon_buffer, lon_range[2] + lon_buffer),
      ylim = c(lat_range[1] - lat_buffer, lat_range[2] + lat_buffer),
      expand = FALSE
    ) +
    labs(
      title = title_text,
      subtitle = paste0(
        region_text, " | Calibrated ", ifelse(var == "abund", "abundance", "biomass"),
        ", ", year_text, "\n", max_obs_text(pdat$density)
      ),
      x = expression(Longitude ~ degree * W),
      y = expression(Latitude ~ degree * N),
      caption = "Calibrated composite: nSPC + belt transect + towed diver. BLT and TOW observations standardized to the nSPC reference method; mean values within 0.1° cells."
    ) +
    theme_pub +
    theme(
      legend.position = c(0.12, 0.15),                            # Coordinates relative to panel (0 to 1)
      legend.justification = c(0, 0),                            # Anchors legend by its bottom-left corner
      legend.background = element_rect(fill = alpha("white", 0.8), color = NA), # Semi-transparent white backing
      legend.key = element_blank(),
      legend.margin = margin(4, 6, 4, 6)
    )
  
  dx <- diff(lon_range) + 2 * lon_buffer
  dy <- diff(lat_range) + 2 * lat_buffer
  map_ratio <- dx / max(dy, 0.1)
  
  # Fixed overhead estimates in inches for text, titles, margins, and axis labels
  v_overhead <- 2.2  # Title, subtitle, caption, x-axis label
  h_overhead <- 1.5  # Y-axis label, tick text, panel margins
  
  # Target slide limits (16x9 Widescreen)
  max_w <- 16.0
  max_h <- 9.0
  
  # Calculate exact panel-constrained dimensions
  if (map_ratio >= (max_w - h_overhead) / (max_h - v_overhead)) {
    # Wide aspect ratio: Constrained by max width (16 in)
    map_b_width  <- max_w
    panel_w      <- max_w - h_overhead
    panel_h      <- panel_w / map_ratio
    map_b_height <- min(max_h, max(4.5, panel_h + v_overhead))
  } else {
    # Tall / Square aspect ratio: Constrained by max height (9 in)
    map_b_height <- max_h
    panel_h      <- max_h - v_overhead
    panel_w      <- panel_h * map_ratio
    map_b_width  <- min(max_w, max(7.0, panel_w + h_overhead))
  }
  
  ggsave(
    filename = paste0("output/plot/map_b_", species[s], "_", var, "_", region_tag, ".png"),
    plot = p, 
    height = map_b_height, 
    width = map_b_width,
    units = "in", 
    dpi = 300, 
    bg = "white"
  )
  
  # -----------------------------------------------------------------------
  # MAP C: composite by island
  # -----------------------------------------------------------------------
  
  pdat <- df %>%
    filter(density > 0, method == "nSPC_BLT_TOW") %>%
    group_by(method, island, longitude, latitude) %>%
    summarise(density = mean(density), .groups = "drop")
  
  pm <- prep_map_bins(pdat)
  pdat <- pm$dat
  n_islands <- length(unique(pdat$island))
  ncol_c <- min(3, n_islands)
  
  p = ggplot(pdat, aes(longitude, latitude)) +
    geom_point(aes(size = density_bin, fill = density_bin),
               shape = 21, alpha = 0.85, colour = "grey20", stroke = 0.25) +
    scale_fill_manual(values = pm$fill_vals, breaks = pm$levels, name = unit, drop = FALSE) +
    scale_size_manual(values = pm$size_vals, breaks = pm$levels, name = unit, drop = FALSE) +
    guides(fill = guide_legend(title.position = "top",
                               override.aes = list(shape = 21, size = unname(pm$size_vals),
                                                   colour = "grey20", alpha = 0.9)),
           size = "none") +
    facet_wrap(~island, scales = "free", ncol = ncol_c) +
    scale_x_continuous(n.breaks = 3, labels = label_number(accuracy = 0.1)) +
    scale_y_continuous(n.breaks = 3, labels = label_number(accuracy = 0.1)) +
    labs(
      title = title_text,
      subtitle = paste0(region_text, " | Calibrated ", ifelse(var == "abund", "abundance", "biomass"),
                        " by island, ", year_text, "\n", max_obs_text(pdat$density)),
      x = expression(Longitude ~ degree * W),
      y = expression(Latitude ~ degree * N),
      caption = "Composite nSPC-equivalent observations. Point size and color represent mean values at each survey location."
    ) +
    theme_pub +
    theme(legend.position = "right")
  
  dims <- facet_dims(n_islands, ncol_c, 3.6, 3.1)
  
  ggsave(
    filename = paste0("output/plot/map_c_", species[s], "_", var, "_", region_tag, ".png"),
    plot = p, height = dims["height"], width = dims["width"],
    units = "in", dpi = 300, bg = "white"
  )
  
  # -----------------------------------------------------------------------
  # MAP D: spatial distribution by survey method and year
  # -----------------------------------------------------------------------
  
  pdat <- df %>%
    filter(density > 0) %>%
    mutate(longitude = round(longitude, 1), latitude = round(latitude, 1)) %>%
    group_by(method, longitude, latitude, year) %>%
    summarise(density = mean(density), .groups = "drop")
  
  pm <- prep_map_bins(pdat)
  pdat <- pm$dat
  n_years <- length(unique(pdat$year))
  n_methods <- length(unique(pdat$method))
  
  p = ggplot(pdat, aes(longitude, latitude)) +
    geom_point(aes(size = density_bin, fill = density_bin),
               shape = 21, alpha = 0.80, colour = "grey20", stroke = 0.2) +
    scale_fill_manual(values = pm$fill_vals, breaks = pm$levels, name = unit, drop = FALSE) +
    scale_size_manual(values = pm$size_vals, breaks = pm$levels, name = unit, drop = FALSE) +
    guides(fill = guide_legend(title.position = "top",
                               override.aes = list(shape = 21, size = unname(pm$size_vals),
                                                   colour = "grey20", alpha = 0.9)),
           size = "none") +
    facet_grid(method ~ year, labeller = labeller(method = method_labs)) +
    labs(
      title = title_text,
      subtitle = paste0(region_text, " | Spatial distribution by survey method and year, ", year_text,
                        "\n", max_obs_text(pdat$density)),
      x = NULL, y = NULL,
      caption = "Point size and color represent mean values within 0.1° cells."
    ) +
    theme_pub +
    theme(
      legend.position = "bottom",
      axis.ticks = element_blank(),
      axis.text = element_blank()
    )
  
  ggsave(
    file = paste0("output/plot/map_d_", species[s], "_", var, "_", region_tag, ".png"),
    plot = p,
    height = max(6, 2.0 + 1.65 * n_methods),
    width = max(10, 2.0 + 1.55 * n_years),
    units = "in", dpi = 300, bg = "white"
  )
  
  # -----------------------------------------------------------------------
  # DEPTH
  # -----------------------------------------------------------------------
  
  p = df %>%
    filter(method == "nSPC_BLT_TOW") %>%
    mutate(depth = round(depth, 1)) %>%
    group_by(method, region, depth) %>%
    summarise(density = mean(density, na.rm = TRUE), .groups = "drop") %>%
    ggplot(aes(depth, density)) +
    geom_point(aes(fill = density), shape = 21, alpha = 0.85,
               size = 3, show.legend = FALSE, colour = "grey20") +
    scale_fill_gradientn(colours = matlab.like(100), trans = "sqrt") +
    labs(
      title = title_text,
      subtitle = paste0(region_text, " | Calibrated ", ifelse(var == "abund", "abundance", "biomass"),
                        " by depth, ", year_text),
      x = "Depth (m)", y = unit
    ) +
    facet_grid(~region) +
    theme_pub
  
  # -----------------------------------------------------------------------
  # TS A: survey-specific temporal series
  # -----------------------------------------------------------------------
  
  pdat <- df %>%
    filter(method != "nSPC_BLT_TOW") %>%
    group_by(year, method, region) %>%
    summarise(mean_density = mean(density),
              se_density = sd(density)/sqrt(n()),
              .groups = "drop")
  
  n_regions <- length(unique(pdat$region))
  ncol_ts <- min(2, n_regions)
  
  p = pdat %>%
    ggplot(aes(x = year, y = mean_density, color = method, group = method)) +
    geom_errorbar(aes(ymin = pmax(0, mean_density - se_density),
                      ymax = mean_density + se_density),
                  width = 0, show.legend = FALSE) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2.5) +
    scale_color_discrete(labels = method_labs[c("nSPC", "BLT", "TOW")], name = NULL) +
    scale_x_continuous(breaks = sort(unique(pdat$year))) +
    labs(
      title = title_text,
      subtitle = paste0(region_text, " | Survey-specific calibrated ",
                        ifelse(var == "abund", "abundance", "biomass"), ", ", year_text),
      x = NULL, y = unit
    ) +
    facet_wrap(~region, scales = "free_y", ncol = ncol_ts,
               labeller = labeller(region = region_labs)) +
    theme_pub +
    theme(
      legend.position = "top",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  dims <- facet_dims(n_regions, ncol_ts, 4.8, 3.5)
  
  ggsave(
    p, file = paste0("output/plot/ts_a_", species[s], "_", var, "_", region_tag, ".png"),
    height = dims["height"], width = dims["width"],
    units = "in", dpi = 300, bg = "white"
  )
  
  # -----------------------------------------------------------------------
  # TS B: pooled temporal series
  # -----------------------------------------------------------------------
  
  pdat <- df %>%
    filter(method != "nSPC_BLT_TOW") %>%
    group_by(year, region) %>%
    summarise(mean_density = mean(density),
              se_density = sd(density)/sqrt(n()),
              .groups = "drop")
  
  p = pdat %>%
    ggplot(aes(x = year, y = mean_density)) +
    geom_errorbar(aes(ymin = pmax(0, mean_density - se_density),
                      ymax = mean_density + se_density),
                  width = 0, colour = "grey35") +
    geom_line(linewidth = 0.7, colour = "grey45") +
    geom_point(size = 3, shape = 21, fill = "white", colour = "grey20") +
    scale_x_continuous(breaks = sort(unique(pdat$year))) +
    labs(
      title = title_text,
      subtitle = paste0(region_text, " | Mean calibrated ",
                        ifelse(var == "abund", "abundance", "biomass"), ", ", year_text),
      x = NULL, y = unit,
      caption = "Annual mean across available calibrated survey observations."
    ) +
    facet_wrap(~region, scales = "free_y", ncol = ncol_ts,
               labeller = labeller(region = region_labs)) +
    theme_pub +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  dims <- facet_dims(n_regions, ncol_ts, 4.8, 3.5)
  
  ggsave(
    p, file = paste0("output/plot/ts_b_", species[s], "_", var, "_", region_tag, ".png"),
    height = dims["height"], width = dims["width"],
    units = "in", dpi = 300, bg = "white"
  )
  
  # -----------------------------------------------------------------------
  # TS C: island-level temporal series
  # -----------------------------------------------------------------------
  
  dfi = df %>%
    filter(method != "nSPC_BLT_TOW") %>%
    group_by(year, island) %>%
    summarise(mean_density = mean(density),
              se_density = sd(density)/sqrt(n()),
              .groups = "drop") %>%
    na.omit()
  
  n_islands <- length(unique(dfi$island))
  ncol_ts_c <- min(3, n_islands)
  
  p = dfi %>%
    ggplot(aes(x = year, y = mean_density)) +
    geom_errorbar(aes(ymin = pmax(0, mean_density - se_density),
                      ymax = mean_density + se_density),
                  width = 0, colour = "grey35") +
    geom_line(linewidth = 0.65, colour = "grey45") +
    geom_point(size = 2.8, shape = 21, fill = "white", colour = "grey20") +
    scale_x_continuous(breaks = sort(unique(dfi$year))) +
    labs(
      title = title_text,
      subtitle = paste0(region_text, " | Island-level calibrated ",
                        ifelse(var == "abund", "abundance", "biomass"), ", ", year_text),
      x = NULL, y = unit
    ) +
    facet_wrap(~island, scales = "free_y", ncol = ncol_ts_c) +
    theme_pub +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  dims <- facet_dims(n_islands, ncol_ts_c, 3.8, 3.0)
  
  ggsave(
    p, file = paste0("output/plot/ts_c_", species[s], "_", var, "_", region_tag, ".png"),
    height = dims["height"], width = dims["width"],
    units = "in", dpi = 300, bg = "white"
  )
}
