read_nspc_summary <- function(obs_set, region = NULL, columns = NULL) {
  if (!obs_set %in% c("I", "IUN")) {
    stop("Supported nSPC observation sets are I and IUN.", call. = FALSE)
  }

  path <- file.path(
    "data", "tanaka_2026",
    paste0("site_summary_SPECIES_", obs_set, ".parquet")
  )
  if (!file.exists(path)) {
    stop("Missing nSPC summary input: ", path, call. = FALSE)
  }

  dataset <- arrow::open_dataset(path, format = "parquet")
  available <- dataset$schema$names
  required <- c(
    "SITEVISITID", "METHOD", "REGION", "ISLAND", "SEC_NAME",
    "REEF_ZONE", "DEPTH_BIN", "OBS_YEAR", "DATE_", "LATITUDE",
    "LONGITUDE", "FISH_GROUP", "PRESENT", "DENSITY_M2",
    "BIOMASS_G_M2", "MEAN_DEPTH_M"
  )
  missing <- setdiff(required, available)
  if (length(missing)) {
    stop("Missing columns in ", path, ": ", paste(missing, collapse = ", "),
         call. = FALSE)
  }

  requested_columns <- if (is.null(columns)) required else columns
  missing <- setdiff(requested_columns, available)
  if (length(missing)) {
    stop("Missing requested columns in ", path, ": ",
         paste(missing, collapse = ", "), call. = FALSE)
  }

  read_columns <- unique(c(requested_columns, "OBS_YEAR", "DATE_"))
  query <- dataset
  if (!is.null(region)) {
    query <- dplyr::filter(query, REGION == !!region)
  }
  query <- dplyr::select(query, dplyr::all_of(read_columns))
  result <- dplyr::collect(query)

  if (!nrow(result)) {
    stop("No nSPC summary rows found for ", obs_set,
         if (!is.null(region)) paste0(" in ", region), call. = FALSE)
  }

  required_values <- intersect(
    c("SITEVISITID", "OBS_YEAR", "FISH_GROUP", "METHOD", "REGION",
      "ISLAND", "REEF_ZONE", "DEPTH_BIN", "DATE_", "LATITUDE",
      "LONGITUDE", "MEAN_DEPTH_M", "PRESENT", "DENSITY_M2",
      "BIOMASS_G_M2"),
    names(result)
  )
  for (column in required_values) {
    values <- result[[column]]
    if (anyNA(values)) {
      stop("Missing values in ", column, " in ", path, call. = FALSE)
    }
    if (column %in% c("SITEVISITID", "OBS_YEAR", "LATITUDE", "LONGITUDE",
                      "MEAN_DEPTH_M", "DENSITY_M2", "BIOMASS_G_M2") &&
        any(!is.finite(values))) {
      stop("Non-finite values in ", column, " in ", path, call. = FALSE)
    }
    if (is.character(values) && any(!nzchar(trimws(values)))) {
      stop("Empty values in ", column, " in ", path, call. = FALSE)
    }
  }

  date_parts <- as.POSIXlt(result$DATE_, tz = "UTC")
  date_year <- date_parts$year + 1900L
  mismatch <- date_year != result$OBS_YEAR
  if (any(mismatch)) {
    recoverable <- date_year[mismatch] == 26L &
      result$OBS_YEAR[mismatch] == 2026L
    if (any(!recoverable)) {
      stop("DATE_ year does not match OBS_YEAR in ", path, call. = FALSE)
    }
    date_parts$year[mismatch] <- result$OBS_YEAR[mismatch] - 1900L
    result$DATE_ <- as.POSIXct(date_parts, tz = "UTC")
  }

  if ("PRESENT" %in% names(result) &&
      any(!result$PRESENT %in% c(0, 1))) {
    stop("Invalid presence values in ", path, call. = FALSE)
  }
  for (column in intersect(c("DENSITY_M2", "BIOMASS_G_M2"), names(result))) {
    if (any(result[[column]] < 0)) {
      stop("Negative values in ", column, " in ", path, call. = FALSE)
    }
  }

  internal_columns <- setdiff(c("OBS_YEAR", "DATE_"), requested_columns)
  result[internal_columns] <- NULL
  result
}
