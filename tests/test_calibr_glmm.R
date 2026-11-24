suppressPackageStartupMessages(library(Calibr))

set.seed(117)
data <- expand.grid(BLOCK = paste0("b", 1:8), GROUP = c("A", "B"),
                    METHOD = c("nSPC", "BLT"), REP = 1:8,
                    stringsAsFactors = FALSE)
data$PRESENCE <- rbinom(nrow(data), 1, ifelse(data$GROUP == "A", 0.65, 0.45))
data$DENSITY <- data$PRESENCE * exp(rnorm(nrow(data), mean = 0, sd = 0.5))
rare <- data[1:4, ]
rare$GROUP <- "RARE"
rare$PRESENCE <- 1
rare$DENSITY <- 1
data <- rbind(data, rare)

result <- run_calibr(data, std_method = "nSPC", stat_model = "GLMM", n_sample = 2)
summary <- result$SUMMARY
for (group in c("A", "B")) {
  rows <- summary[summary$GROUP == group, , drop = FALSE]
  stopifnot(nrow(rows) == 2L, all(rows$GCF.SOURCE == "Group-specific"),
            all(is.finite(rows$PRES.GCF)), all(is.finite(rows$POS.GCF)))
}
rare_rows <- summary[summary$GROUP == "RARE", , drop = FALSE]
stopifnot(nrow(rare_rows) > 0L, all(rare_rows$GCF.SOURCE == "Global"))
cat("Calibr GLMM group-specific and global fallback checks passed.\n")
