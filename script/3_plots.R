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
library(patchwork)

rm(list = ls())

source("script/nspc_summary.R")

plot_args <- commandArgs(trailingOnly = TRUE)
var <- if (length(plot_args) >= 1L) plot_args[[1]] else "abund"
region <- if (length(plot_args) >= 2L) plot_args[[2]] else "SAMOA"
plot_selection <- if (length(plot_args) >= 3L) plot_args[[3]] else "all"

if (length(plot_args) > 4L) {
  stop("Usage: 3_plots.R [METRIC [REGION [all|maps [SPECIES]]]]", call. = FALSE)
}
if (!var %in% c("abund", "biom")) {
  stop("Metric must be 'abund' or 'biom'.", call. = FALSE)
}
if (!region %in% c("MHI", "MARIAN", "NWHI", "PRIAs", "SAMOA")) {
  stop("Unsupported region: ", region, call. = FALSE)
}
if (!plot_selection %in% c("all", "maps")) {
  stop("Unsupported plot selection: ", plot_selection, call. = FALSE)
}

available_species = c(
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
)
species <- if (length(plot_args) >= 4L) plot_args[[4]] else "APVI"
if (!species %in% available_species) {
  stop("Unsupported species: ", species, call. = FALSE)
}

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

theme_pub <- theme_minimal(base_size = 12) +
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
  brks <- sort(unique(c(0, pmin(brks, cap), cap)))
  if (length(brks) < 3) brks <- sort(unique(c(0, pretty(c(0, cap), n = 4), cap)))
  accuracy <- 10^floor(log10(min(diff(brks))))

  bin_labs <- paste0(
    scales::number(brks[-length(brks)], accuracy = accuracy), "–",
    scales::number(brks[-1], accuracy = accuracy)
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

longitude_labels <- function(x) {
  values <- ifelse(x > 180, 360 - x, x)
  paste0(scales::label_number()(values), ifelse(x > 180, "°W", "°E"))
}

coast_path <- "data/ne_10m_land.zip"
if (!file.exists(coast_path)) {
  stop("Missing Natural Earth coastline data: ", coast_path, call. = FALSE)
}
coast <- sf::st_geometry(sf::st_read(
  paste0("/vsizip/", normalizePath(coast_path, winslash = "/", mustWork = TRUE),
         "/ne_10m_land.shp"), quiet = TRUE
))
coast_for_extent <- function(xlim, ylim) {
  previous_s2 <- sf::sf_use_s2()
  on.exit(suppressMessages(sf::sf_use_s2(previous_s2)))
  suppressMessages(sf::sf_use_s2(FALSE))
  intervals <- list(c(max(0, xlim[1]), min(180, xlim[2])),
                    c(max(180, xlim[1]), min(360, xlim[2])))
  pieces <- lapply(seq_along(intervals), function(i) {
    interval <- intervals[[i]]
    if (interval[1] >= interval[2]) return(NULL)
    if (i == 2L) interval <- interval - 360
    bbox <- sf::st_bbox(c(xmin = interval[1], ymin = ylim[1],
                          xmax = interval[2], ymax = ylim[2]),
                        crs = sf::st_crs(4326))
    clipped <- suppressMessages(sf::st_crop(coast, bbox))
    if (!length(clipped)) return(NULL)
    polygons <- if (inherits(clipped, c("sfc_GEOMETRY", "sfc_GEOMETRYCOLLECTION"))) {
      sf::st_collection_extract(clipped, "POLYGON")
    } else {
      clipped
    }
    if (!length(polygons)) return(NULL)
    polygons <- sf::st_cast(polygons, "MULTIPOLYGON")
    coords <- as.data.frame(sf::st_coordinates(polygons))
    coords %>% transmute(
      longitude = ifelse(X < 0, X + 360, X), latitude = Y,
      group = paste(i, L3, L2, sep = "-"),
      subgroup = paste(i, L3, L2, L1, sep = "-")
    )
  })
  if (all(vapply(pieces, is.null, logical(1)))) {
    return(data.frame(longitude = numeric(), latitude = numeric(),
                      group = character(), subgroup = character()))
  }
  bind_rows(pieces)
}

island_satellite <- function(xlim, ylim) {
  if (length(xlim) != 2 || length(ylim) != 2 ||
      any(!is.finite(c(xlim, ylim))) || xlim[1] >= xlim[2] || ylim[1] >= ylim[2] ||
      xlim[1] < 0 || xlim[2] > 360 || ylim[1] < -90 || ylim[2] > 90) {
    stop("Island map extent is outside satellite imagery coverage.", call. = FALSE)
  }
  cols <- seq(floor(xlim[1]), ceiling(xlim[2]) - 1)
  rows <- seq(floor(ylim[1]), ceiling(ylim[2]) - 1)
  if (length(cols) * length(rows) > 64) {
    stop("Island map extent requires more than 64 satellite imagery tiles.", call. = FALSE)
  }
  cache_dir <- file.path("data", ".eox_2025_cache")
  if (!dir.exists(cache_dir) && !dir.create(cache_dir, recursive = TRUE)) {
    stop("Cannot create satellite imagery cache: ", cache_dir, call. = FALSE)
  }
  read_tile <- function(path) {
    signature <- readBin(path, "raw", n = 8)
    if (identical(signature[1:3], as.raw(c(0xff, 0xd8, 0xff)))) {
      image <- jpeg::readJPEG(path, native = TRUE)
    } else if (identical(signature, as.raw(c(0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a)))) {
      image <- png::readPNG(path, native = TRUE)
    } else {
      stop("Invalid satellite imagery response: ", path, call. = FALSE)
    }
    if (length(dim(image)) != 2L || nrow(image) != ncol(image) ||
        !nrow(image) %in% c(256L, 512L, 1024L, 2048L)) {
      stop("Invalid satellite imagery tile dimensions: ", path, call. = FALSE)
    }
    image
  }
  load_tile <- function(lon, lat) {
    path <- file.path(cache_dir, sprintf("lon%d_lat%d.tile", lon, lat))
    if (!file.exists(path)) {
      signed_lon <- if (lon >= 180) lon - 360 else lon
      tmp <- tempfile(tmpdir = cache_dir, fileext = ".jpg")
      on.exit(unlink(tmp), add = TRUE)
      for (resolution in c(2048L, 1024L, 512L, 256L)) {
        url <- sprintf(
          paste0("https://tiles.maps.eox.at/wms?service=WMS&version=1.1.1&request=GetMap",
                 "&layers=s2cloudless-2025&styles=&format=image%%2Fjpeg&srs=EPSG%%3A4326",
                 "&bbox=%d,%d,%d,%d&width=%d&height=%d"),
          signed_lon, lat, signed_lon + 1, lat + 1, resolution, resolution
        )
        response <- curl::curl_fetch_disk(url, tmp)
        if (response$status_code == 200L) break
        no_tiles <- response$status_code == 404L &&
          grepl("x-mapcache-error: no tiles containing image data could be retrieved to create map",
                tolower(rawToChar(response$headers)), fixed = TRUE)
        if (!no_tiles || resolution == 256L) {
          stop("Could not download satellite imagery tile ", lon, ", ", lat,
               " (HTTP ", response$status_code, ").", call. = FALSE)
        }
      }
      image <- read_tile(tmp)
      if (nrow(image) != resolution) {
        stop("Unexpected satellite imagery tile dimensions: ", tmp, call. = FALSE)
      }
      if (!file.rename(tmp, path)) stop("Cannot save satellite imagery tile: ", path, call. = FALSE)
    } else {
      image <- read_tile(path)
    }
    resolution <- nrow(image)
    if (resolution < 2048L) {
      warning(sprintf("Satellite imagery is lower resolution in cell %d, %d.", lon, lat),
              call. = FALSE)
    }
    col_min <- max(0, floor((xlim[1] - lon) * resolution))
    col_max <- min(resolution, ceiling((xlim[2] - lon) * resolution))
    row_min <- max(0, floor((lat + 1 - ylim[2]) * resolution))
    row_max <- min(resolution, ceiling((lat + 1 - ylim[1]) * resolution))
    # Native rasters are stored row-major, unlike R matrices.
    cropped <- image[seq.int(col_min + 1, col_max),
                     seq.int(row_min + 1, row_max), drop = FALSE]
    dim(cropped) <- c(row_max - row_min, col_max - col_min)
    if (attr(image, "channels") == 4L && any(bitwShiftR(cropped, 24L) == 0L)) {
      warning(sprintf("Satellite imagery has missing coverage in cell %d, %d.", lon, lat),
              call. = FALSE)
    }
    class(cropped) <- "nativeRaster"
    attr(cropped, "channels") <- attr(image, "channels")
    annotation_raster(
      cropped, xmin = lon + col_min / resolution, xmax = lon + col_max / resolution,
      ymin = lat + 1 - row_max / resolution, ymax = lat + 1 - row_min / resolution,
      interpolate = TRUE
    )
  }
  unlist(lapply(rows, function(lat) lapply(cols, function(lon) load_tile(lon, lat))),
         recursive = FALSE)
}

for (s in 1:length(species)) {

  if (plot_selection == "maps") {
    data_path <- paste0("output/calibr_df/calibr_", species[s], "_", var, "_", region, ".RData")
    if (!file.exists(data_path)) {
      stop("Missing calibrated data for maps: ", data_path, call. = FALSE)
    }
    data_env <- new.env(parent = emptyenv())
    load(data_path, envir = data_env)
    if (!exists("df", envir = data_env, inherits = FALSE) ||
        !all(c("region", "island", "method", "longitude", "latitude", "density", "year") %in% names(data_env$df))) {
      stop("Invalid calibrated data for maps: ", data_path, call. = FALSE)
    }
    df <- data_env$df
  } else {
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
    
    spc_density <- if (var == "abund") "DENSITY_M2" else "BIOMASS_G_M2"
    spc_data <- read_nspc_summary(
      "I", region[r],
      columns = c("ISLAND", "METHOD", "DATE_", "LATITUDE", "LONGITUDE",
                  "FISH_GROUP", "MEAN_DEPTH_M", "PRESENT", spc_density)
    )
    spc <- spc_data %>%
      filter(FISH_GROUP == species[s]) %>%
      transmute(
        ISLAND, DEPTH = MEAN_DEPTH_M, METHOD, DATE_, LATITUDE, LONGITUDE,
        SPECIES = FISH_GROUP, DENSITY = .data[[spc_density]],
        PRESENCE = as.integer(PRESENT)
      )
    
    if (dim(spc)[1] == 0) {
      spc <- spc_data %>%
        distinct(ISLAND, DEPTH = MEAN_DEPTH_M, METHOD, DATE_, LATITUDE, LONGITUDE) %>%
        mutate(SPECIES = species[s], DENSITY = 0, PRESENCE = 0)
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
  }
  
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
  wide_region <- diff(lon_range) > 2 * diff(lat_range) && n_methods > 1
  map_cols <- if (wide_region) 1 else n_methods
  region_coast <- coast_for_extent(
    c(lon_range[1] - lon_buffer, lon_range[2] + lon_buffer),
    c(lat_range[1] - lat_buffer, lat_range[2] + lat_buffer)
  )

  p = ggplot(pdat, aes(longitude, latitude)) +
    geom_polygon(
      data = region_coast,
      aes(x = longitude, y = latitude, group = group, subgroup = subgroup),
      fill = "grey85", colour = "grey60", linewidth = 0.25,
      rule = "evenodd", inherit.aes = FALSE
    ) +
    geom_point(aes(size = density_bin, fill = density_bin),
               shape = 21, alpha = 0.85, colour = "grey20", stroke = 0.25) +
    scale_fill_manual(values = pm$fill_vals, breaks = pm$levels, name = unit, drop = FALSE) +
    scale_size_manual(values = pm$size_vals, breaks = pm$levels, name = unit, drop = FALSE) +
    guides(fill = guide_legend(title.position = "top",
                               override.aes = list(shape = 21, size = unname(pm$size_vals),
                                                   colour = "grey20", alpha = 0.9)),
           size = "none") +
    facet_wrap(~method, ncol = map_cols, labeller = labeller(method = method_labs)) +
    scale_x_continuous(n.breaks = 3, labels = longitude_labels) +
    coord_equal(xlim = c(lon_range[1] - lon_buffer, lon_range[2] + lon_buffer),
                ylim = c(lat_range[1] - lat_buffer, lat_range[2] + lat_buffer),
                expand = FALSE) +
    labs(
      title = title_text,
      subtitle = paste0(region_text, " | Survey-specific calibrated ", ifelse(var == "abund", "abundance", "biomass"),
                        ", ", year_text, "\n", max_obs_text(pdat$density)),
      x = "Longitude",
      y = "Latitude (degrees)",
      caption = paste("BLT and TOW observations standardized to the nSPC reference method.",
                      "Mean values within 0.1° cells; values above the 99th percentile capped for display.",
                      sep = "\n")
    ) +
    theme_pub +
    theme(legend.position = "right", panel.spacing.x = grid::unit(1.5, "lines"),
          panel.background = element_rect(fill = "#edf5f7", colour = NA))

  map_ratio <- diff(lon_range) / max(diff(lat_range), 0.1)
  dims <- if (wide_region) {
    c(width = 11, height = 9)
  } else if (map_ratio < 1) {
    c(width = max(8, min(15, 3.5 + n_methods * 4 * map_ratio)), height = 8)
  } else {
    facet_dims(n_methods, n_methods, 4.3, 4.1)
  }
  
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
  region_coast <- coast_for_extent(
    c(lon_range[1] - lon_buffer, lon_range[2] + lon_buffer),
    c(lat_range[1] - lat_buffer, lat_range[2] + lat_buffer)
  )
  
  p <- ggplot() +
    geom_polygon(
      data = region_coast,
      aes(x = longitude, y = latitude, group = group, subgroup = subgroup),
      fill = "grey85", colour = "grey60", linewidth = 0.25, rule = "evenodd"
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
    scale_x_continuous(labels = longitude_labels) +
    coord_equal(
      xlim = c(lon_range[1] - lon_buffer, lon_range[2] + lon_buffer),
      ylim = c(lat_range[1] - lat_buffer, lat_range[2] + lat_buffer),
      expand = FALSE
    ) +
    labs(
      title = title_text,
      subtitle = paste0(
        region_text, " | Regional overview of calibrated ",
        ifelse(var == "abund", "abundance", "biomass"), ", ", year_text
      ),
      x = "Longitude",
      y = "Latitude (degrees)",
      caption = paste("Calibrated composite: nSPC + belt transect + towed diver.",
                      "Points show mean positive observations within 0.1° cells, capped at the 99th percentile.",
                      "The island panels below instead include zero observations and have their own legend.",
                      sep = "\n")
    ) +
    theme_pub +
    theme(
      legend.position = "right",
      legend.key = element_blank(),
      legend.margin = margin(4, 6, 4, 6),
      panel.background = element_rect(fill = "#edf5f7", colour = NA)
    )
  
  dx <- diff(lon_range) + 2 * lon_buffer
  dy <- diff(lat_range) + 2 * lat_buffer
  map_ratio <- dx / max(dy, 0.1)
  
  # Fixed overhead estimates in inches for text, titles, margins, and axis labels
  v_overhead <- 2.2  # Title, subtitle, caption, x-axis label
  h_overhead <- 3.5  # Y-axis label, tick text, panel margins, and legend
  
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
  
  regional_map <- p
  
  # -----------------------------------------------------------------------
  # MAP C: composite by island
  # -----------------------------------------------------------------------
  
  pdat <- df %>%
    filter(method == "nSPC_BLT_TOW") %>%
    group_by(island, longitude, latitude) %>%
    summarise(density = mean(density), .groups = "drop")
  if (!nrow(pdat) || any(!is.finite(pdat$density)) || !any(pdat$density > 0)) {
    stop("Island map requires finite composite observations and positive values.",
         call. = FALSE)
  }
  pm <- prep_map_bins(pdat %>% filter(density > 0))
  pm$size_vals <- pm$size_vals * 0.55
  positive <- pm$dat
  island_names <- sort(unique(pdat$island))
  views <- unlist(lapply(island_names, function(island_name) {
    sites <- pdat %>% filter(island == island_name)
    full_x <- range(sites$longitude)
    full_y <- range(sites$latitude)
    result <- list(list(island = island_name, x = full_x, y = full_y, zoom = FALSE))
    if (nrow(sites) >= 20) {
      zoom_x <- unname(quantile(sites$longitude, c(0.05, 0.95)))
      zoom_y <- unname(quantile(sites$latitude, c(0.05, 0.95)))
      if ((diff(zoom_x) > 0 && diff(full_x) > 3 * diff(zoom_x)) ||
          (diff(zoom_y) > 0 && diff(full_y) > 3 * diff(zoom_y))) {
        result <- c(result, list(list(island = island_name, x = zoom_x,
                                     y = zoom_y, zoom = TRUE)))
      }
    }
    result
  }), recursive = FALSE)
  ncol_c <- min(3, length(views))

  island_plots <- lapply(views, function(island_view) {
    island_name <- island_view$island
    island_sites <- pdat %>% filter(island == island_name)
    island_positive <- positive %>% filter(island == island_name)
    xr <- island_view$x
    yr <- island_view$y
    x_pad <- max(diff(xr) * 0.12, 0.002)
    y_pad <- max(diff(yr) * 0.12, 0.002)
    xlim <- xr + c(-x_pad, x_pad)
    ylim <- yr + c(-y_pad, y_pad)
    visible <- sum(island_sites$longitude >= xlim[1] & island_sites$longitude <= xlim[2] &
                   island_sites$latitude >= ylim[1] & island_sites$latitude <= ylim[2])

    p_island <- ggplot() + island_satellite(xlim, ylim)
    p_island +
      geom_point(data = island_sites %>% filter(density == 0),
                 aes(longitude, latitude), shape = 21, fill = "white",
                 colour = "grey20", size = 1.1, stroke = 0.25, alpha = 0.8) +
      geom_point(data = island_positive,
                 aes(longitude, latitude, size = density_bin, fill = density_bin),
                 shape = 21, alpha = 0.9, colour = "white", stroke = 0.45) +
      scale_fill_manual(values = pm$fill_vals, breaks = pm$levels, name = unit, drop = FALSE) +
      scale_size_manual(values = pm$size_vals, breaks = pm$levels, name = unit, drop = FALSE) +
      guides(fill = "none", size = "none") +
      scale_x_continuous(n.breaks = 3, labels = longitude_labels,
                         guide = guide_axis(check.overlap = TRUE)) +
      scale_y_continuous(n.breaks = 3, labels = label_number(),
                         guide = guide_axis(check.overlap = TRUE)) +
      coord_quickmap(xlim = xlim, ylim = ylim, expand = FALSE) +
      labs(title = paste0(island_name, if (island_view$zoom) " (close-up)" else "",
                          "  |  ", visible, " locations"),
           x = NULL, y = NULL) +
      theme_minimal(base_size = 10) +
      theme(
        plot.title = element_text(face = "bold", size = 10, hjust = 0.5,
                                  margin = margin(b = 5)),
        panel.background = element_rect(fill = "#edf5f7", colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "white", linewidth = 0.3),
        axis.text = element_text(size = 8, colour = "grey30")
      )
  })

  legend_plot <- ggplot(positive, aes(longitude, latitude)) +
    geom_point(aes(size = density_bin, fill = density_bin),
               shape = 21, colour = "grey20", alpha = 0.9, show.legend = TRUE) +
    scale_fill_manual(values = pm$fill_vals, breaks = pm$levels, name = unit, drop = FALSE) +
    scale_size_manual(values = pm$size_vals, breaks = pm$levels, name = unit, drop = FALSE) +
    guides(fill = guide_legend(title.position = "top",
                               override.aes = list(shape = 21, size = unname(pm$size_vals),
                                                   colour = "grey20", alpha = 0.9)),
           size = "none") +
    theme_minimal(base_size = 10) +
    theme(legend.position = "right")
  legend_table <- ggplotGrob(legend_plot)
  legend_index <- which(legend_table$layout$name == "guide-box-right")
  if (length(legend_index) != 1L) stop("Could not build island map legend.", call. = FALSE)
  legend <- legend_table$grobs[[legend_index]]
  island_grobs <- lapply(island_plots, function(island_plot) {
    wrap_elements(full = ggplotGrob(island_plot))
  })

  p <- (wrap_plots(island_grobs, ncol = ncol_c) | wrap_elements(full = legend)) +
    plot_layout(widths = c(5, 1)) +
    plot_annotation(
      title = title_text,
      subtitle = paste0(region_text, " | Mean calibrated ",
                        ifelse(var == "abund", "abundance", "biomass"),
                        " at surveyed locations, ", year_text),
      caption = paste("Each point is the mean of all composite observations at that location, including zeros.",
                      "White points are zero means; color and size show positive means, capped at the 99th percentile.",
                      "Panels have independent geographic scales; close-ups supplement full island views.",
                      "EOxCloudless https://cloudless.eox.at by EOX IT Services GmbH https://eox.at",
                      "(Contains modified Copernicus Sentinel data 2025); imagery resampled, points added.",
                      "CC BY-NC-SA 4.0: https://creativecommons.org/licenses/by-nc-sa/4.0/",
                      "Historical imagery may have gaps and is not navigation-grade.",
                      sep = "\n"),
      theme = theme_pub
    )

  dims <- facet_dims(length(views), ncol_c, 4.1, 3.7)
  
  ggsave(
    filename = paste0("output/plot/map_c_", species[s], "_", var, "_", region_tag, ".png"),
    plot = p, height = dims["height"], width = dims["width"],
    units = "in", dpi = 300, bg = "white"
  )

  island_detail <- p + plot_annotation(
    title = "Island detail",
    subtitle = paste0("Mean calibrated ", ifelse(var == "abund", "abundance", "biomass"),
                      " at surveyed locations, including zeros"),
    caption = paste("Each point is the mean of all composite observations at that location, including zeros.",
                    "White points are zero means; color and size show positive means, capped at the 99th percentile.",
                    "Panels have independent geographic scales; close-ups supplement full island views.",
                    "EOxCloudless https://cloudless.eox.at by EOX IT Services GmbH https://eox.at",
                    "(Contains modified Copernicus Sentinel data 2025); imagery resampled, points added.",
                    "CC BY-NC-SA 4.0: https://creativecommons.org/licenses/by-nc-sa/4.0/",
                    "Historical imagery may have gaps and is not navigation-grade.",
                    sep = "\n"),
    theme = theme_pub
  )
  map_b <- wrap_elements(full = ggplotGrob(regional_map)) /
    wrap_elements(full = patchworkGrob(island_detail)) +
    plot_layout(heights = c(map_b_height, dims["height"]))
  ggsave(
    filename = paste0("output/plot/map_b_", species[s], "_", var, "_", region_tag, ".png"),
    plot = map_b, height = map_b_height + dims["height"],
    width = max(map_b_width, dims["width"]), units = "in", dpi = 300, bg = "white"
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
  unknown_methods <- setdiff(unique(pdat$method), names(method_labs))
  if (length(unknown_methods)) {
    stop("Unknown map D survey methods: ", paste(unknown_methods, collapse = ", "),
         call. = FALSE)
  }
  methods_d <- names(method_labs)[names(method_labs) %in% unique(pdat$method)]
  ncol_d <- 3
  pm$size_vals <- pm$size_vals * 0.6
  lon_range <- range(pdat$longitude)
  lat_range <- range(pdat$latitude)
  xlim_d <- lon_range + c(-1, 1) * max(diff(lon_range) * 0.05, 0.002)
  ylim_d <- lat_range + c(-1, 1) * max(diff(lat_range) * 0.05, 0.002)
  region_coast <- coast_for_extent(xlim_d, ylim_d)
  rows_d <- vapply(methods_d, function(method_name) {
    ceiling(n_distinct(pdat$year[pdat$method == method_name]) / ncol_d)
  }, numeric(1))

  method_plots <- lapply(methods_d, function(method_name) {
    ggplot(pdat %>% filter(method == method_name), aes(longitude, latitude)) +
      geom_polygon(data = region_coast,
                   aes(longitude, latitude, group = group, subgroup = subgroup),
                   fill = "grey85", colour = "grey60", linewidth = 0.2,
                   rule = "evenodd", inherit.aes = FALSE) +
      geom_point(aes(size = density_bin, fill = density_bin),
                 shape = 21, alpha = 0.80, colour = "grey20", stroke = 0.2,
                 show.legend = TRUE) +
      scale_fill_manual(values = pm$fill_vals, breaks = pm$levels, name = unit, drop = FALSE) +
      scale_size_manual(values = pm$size_vals, breaks = pm$levels, name = unit, drop = FALSE) +
      guides(fill = guide_legend(title.position = "top",
                                 override.aes = list(shape = 21, size = unname(pm$size_vals),
                                                     colour = "grey20", alpha = 0.9)),
             size = "none") +
      facet_wrap(vars(year), ncol = ncol_d) +
      scale_x_continuous(n.breaks = 3, labels = longitude_labels,
                         guide = guide_axis(check.overlap = TRUE)) +
      scale_y_continuous(n.breaks = 3, guide = guide_axis(check.overlap = TRUE)) +
      coord_cartesian(xlim = xlim_d, ylim = ylim_d, expand = FALSE) +
      labs(title = unname(method_labs[method_name]), x = NULL, y = NULL) +
      theme_pub +
      theme(
        plot.title = element_text(face = "bold", size = 13, margin = margin(b = 6)),
        plot.margin = margin(6, 8, 8, 8),
        axis.text = element_text(size = 8, colour = "grey30"),
        panel.spacing = grid::unit(0.7, "lines"),
        panel.background = element_rect(fill = "#edf5f7", colour = NA)
      )
  })

  p <- wrap_plots(method_plots, ncol = 1, heights = rows_d + 0.35,
                  guides = "collect") &
    theme(legend.position = "bottom")
  p <- p + plot_annotation(
    title = title_text,
    subtitle = paste0(region_text, " | Spatial distribution by survey method and year, ", year_text,
                      "\n", sub("Maximum observed", "Maximum cell mean",
                                max_obs_text(pdat$density), fixed = TRUE)),
    caption = paste("Each section is a survey method; panels show years with positive observations.",
                    "All panels share the same geographic extent and color/size scale.",
                    "Points are means of positive observations within 0.1° cells, capped at the 99th percentile for display.",
                    "Small or widely separated islands are easier to inspect in the island-detail maps.",
                    sep = "\n"),
    theme = theme_pub + theme(legend.position = "bottom")
  )

  dims <- c(width = 14, height = 3 + sum(rows_d) * 1.75 + length(methods_d) * 0.45)
  ggsave(
    file = paste0("output/plot/map_d_", species[s], "_", var, "_", region_tag, ".png"),
    plot = p,
    height = dims["height"],
    width = dims["width"],
    units = "in", dpi = 300, bg = "white"
  )

  if (plot_selection == "maps") next
  
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
