windows <- read.csv("tables/climwin_recovered/candidate_windows.csv")
randomised <- read.csv("tables/climwin_recovered/randomised_windows.csv")
w <- subset(windows, WindowOpen - WindowClose >= 2)
best <- w[which.min(w$deltaAICc), ]
par(mfrow = c(1,2), mar = c(4,4,2,1))
colours <- hcl.colors(100, "Viridis")
shade <- cut(w$deltaAICc, breaks = 100, labels = FALSE)
plot(w$WindowOpen, w$WindowClose, pch = 15, cex = 1.4, col = colours[shade],
     xlab = "Opening offset (weeks)", ylab = "Closing offset (weeks)", main = "Eligible windows")
points(best$WindowOpen, best$WindowClose, pch = 8, cex = 1.7)
hist(randomised$deltaAICc, xlim = range(best$deltaAICc, randomised$deltaAICc),
     main = "500 saved full-search randomisations", xlab = "Delta AICc", col = "#9CC5B3", border = "white")
abline(v = best$deltaAICc, col = "#D55E00", lwd = 2)
