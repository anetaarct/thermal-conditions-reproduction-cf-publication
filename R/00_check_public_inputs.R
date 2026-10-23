# New public-input check; not the historical download/cleaning script.
# Run from the publication repository root.
library(readr)
library(dplyr)
daily <- read_csv("data_derived/smhi_hoburg_daily_temperature.csv", show_col_types = FALSE)
annual <- read_csv("data_derived/annual_main_window_temperature.csv", show_col_types = FALSE)
phenology <- read_csv("data_derived/annual_raw_median_laying_date.csv", show_col_types = FALSE)
stopifnot(nrow(daily) == 16802L, !anyDuplicated(daily$date))
stopifnot(nrow(annual) == 46L, all(annual$n_days == 43L))
stopifnot(nrow(phenology) == 45L, !2020L %in% phenology$year)
rebuilt <- daily |>
  filter(format(date, "%m-%d") >= "03-26", format(date, "%m-%d") <= "05-07") |>
  group_by(year) |>
  summarise(temp_mean = mean(temperature), n_days = n(), .groups = "drop")
check <- inner_join(annual, rebuilt, by = "year", suffix = c("_saved", "_rebuilt"))
stopifnot(nrow(check) == 46L, all(check$n_days_rebuilt == 43L))
stopifnot(max(abs(check$temp_mean_saved - check$temp_mean_rebuilt)) < 1e-10)
print(check)
