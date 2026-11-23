# Inflates every cost in inputs/parameters.csv to the cost year in
# inputs/settings.csv (cost_year) with the CPI-U medical care index
# (data-raw/raw/cpi_medical_care.csv; BLS CUUS0000SAM annual averages).
# Rows with price_value and price_year are rewritten: value = price_value x
# CPI[cost_year] / CPI[price_year]; se, low and high are scaled by the same factor.
# Run from the project folder:  Rscript data-raw/inflate_costs.R
suppressPackageStartupMessages(library(readr))
P <- read_csv("inputs/parameters.csv", col_types = cols(.default = "c"), na = "")
S <- read_csv("inputs/settings.csv", col_types = cols(.default = "c"))
cpi <- read_csv("data-raw/raw/cpi_medical_care.csv", show_col_types = FALSE)
base <- as.integer(S$value[S$setting == "cost_year"]); if (!length(base)) base <- max(cpi$year)
idx <- function(y) { v <- cpi$index[cpi$year == y]; if (!length(v)) stop("no CPI for year ", y); v }
fmt <- function(x) ifelse(is.na(x), NA, trimws(formatC(x, digits = 15, format = "fg")))
rows <- which(!is.na(P$price_value) & !is.na(P$price_year))
for (i in rows) {
  y <- as.integer(P$price_year[i]); f <- idx(min(y, base)) ; f <- idx(base) / f     # years after base: no deflation index, used as is
  v0 <- as.numeric(P$value[i]); v <- round(as.numeric(P$price_value[i]) * f)
  r <- if (is.finite(v0) && v0 > 0) v / v0 else f                                 # rescale spread with the value
  P$value[i] <- fmt(v)
  for (k in c("se", "low", "high")) if (!is.na(P[[k]][i])) P[[k]][i] <- fmt(signif(as.numeric(P[[k]][i]) * r, 6))
}
write_csv(P, "inputs/parameters.csv", na = "")
message(sprintf("inflated %d costs to %d USD", length(rows), base))
