# Run from the repository root after 0_prep_nspc_data.R.
# Site totals sum all stored species, size bins, and depth rows, matching
# the site-level aggregation in 1_run_calibr.R.
library(data.table)
library(ggplot2)

regions <- c("MARIAN", "MHI", "NWHI", "PRIAs", "SAMOA")
versions <- c("Previous (2009-2022)" = "20092022",
              "Updated (2009-2024)" = "20092024")
metrics <- c(abund = "Abundance (fish/m2)", biom = "Biomass (g/m2)")
summaries <- list()

for (region in regions) {
  for (metric in names(metrics)) {
    for (version in names(versions)) {
      path <- sprintf("data/nSPC.site.%s.size.%s.%s.rds",
                      metric, versions[[version]], region)
      if (!file.exists(path)) stop("Missing input: ", path)
      message("Summarizing ", path)
      dat <- as.data.table(readRDS(path))
      value_col <- paste0(metric, ".site")
      required <- c("SITEVISITID", "OBS_YEAR", value_col)
      if (!all(required %in% names(dat))) stop("Missing columns in: ", path)
      if (anyNA(dat$SITEVISITID) || anyNA(dat$OBS_YEAR) ||
          any(!is.finite(dat[[value_col]])) || any(dat[[value_col]] < 0)) {
        stop("Missing site/year or invalid density values in: ", path)
      }
      site_totals <- dat[, .(density = sum(get(value_col))),
                         by = .(SITEVISITID, OBS_YEAR)]
      annual <- site_totals[, .(mean_density = mean(density),
                                n_sites = .N), by = OBS_YEAR]
      annual[, `:=`(region = region, metric = metrics[[metric]],
                    version = version)]
      summaries[[length(summaries) + 1L]] <- annual
      rm(dat, site_totals)
      gc()
    }
  }
}

comparison <- rbindlist(summaries)
comparison[, version := factor(version, levels = names(versions))]
comparison[, region := factor(region, levels = regions)]
comparison[, metric := factor(metric, levels = unname(metrics))]
setorder(comparison, metric, region, version, OBS_YEAR)

plot <- ggplot(comparison,
               aes(OBS_YEAR, mean_density, color = version, linetype = version)) +
  geom_vline(xintercept = 2022.5, color = "grey65", linetype = "dotted") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8, aes(shape = version)) +
  facet_wrap(vars(metric, region), ncol = 5, scales = "free_y") +
  scale_color_manual(values = c("#555555", "#0072B2")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  scale_shape_manual(values = c(1, 16)) +
  scale_x_continuous(breaks = c(2009, 2014, 2019, 2024)) +
  labs(title = "nSPC: previous vs. updated datasets",
       subtitle = "Annual mean of all-species site totals; dotted line marks the newly added years",
       x = "Survey year", y = "Mean site density",
       color = NULL, linetype = NULL, shape = NULL,
       caption = paste("Totals sum stored species, size-bin and depth rows (as in calibration).",
                       "Overlapping years may differ because the raw source was updated.",
                       "Panel y-axis scales differ.")) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank())

if (!dir.exists("output")) dir.create("output")
ggsave("output/nspc_previous_vs_updated.png", plot,
       width = 16, height = 7, dpi = 180)
print(plot)
message("Saved output/nspc_previous_vs_updated.png")
