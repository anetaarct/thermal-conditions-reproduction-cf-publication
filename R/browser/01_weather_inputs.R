daily <- read.csv("data_derived/smhi_hoburg_daily_temperature.csv")
annual <- read.csv("data_derived/annual_main_window_temperature.csv")
laying <- read.csv("data_derived/annual_raw_median_laying_date.csv")
daily$date <- as.Date(daily$date)
stopifnot(nrow(daily) == 16802, !anyDuplicated(daily$date),
          nrow(annual) == 46, nrow(laying) == 45)
keep <- format(daily$date, "%m-%d") >= "03-26" &
        format(daily$date, "%m-%d") <= "05-07"
x <- daily[keep, ]
rebuilt <- aggregate(temperature ~ year, x, mean)
days <- table(x$year)
stopifnot(all(days == 43),
  max(abs(rebuilt$temperature - annual$temp_mean[match(rebuilt$year, annual$year)])) < 1e-10)
print(head(annual))
cat("Daily records:", nrow(daily), "\nYears:", nrow(rebuilt),
    "\nDays per window:", unique(as.integer(days)), "\nInput checks passed.\n")
