## PART 1 — growth_data
## Computes summary stats, boxplots, 10yr growth, t-test

infile  <- "data/growth_data.csv"
outdirT <- "output"
outdirF <- "figures"

growth <- read.csv(infile, stringsAsFactors = FALSE)

# Q6: column names
write.csv(as.data.frame(names(growth)), file.path(outdirT, "q06_column_names.csv"), row.names = FALSE)

# Q7: mean & sd at start (2005) and end (2020) by Site
mean_summary <- aggregate(cbind(Circumf_2005_cm, Circumf_2020_cm) ~ Site, data = growth, FUN = mean)
sd_summary   <- aggregate(cbind(Circumf_2005_cm, Circumf_2020_cm) ~ Site, data = growth, FUN = sd)
summary_table <- merge(mean_summary, sd_summary, by = "Site", suffixes = c("_Mean", "_SD"))
write.csv(summary_table, file.path(outdirT, "q07_mean_sd_2005_2020_by_site.csv"), row.names = FALSE)

# Q8: boxplots (start vs end by Site)
png(file.path(outdirF, "fig_q08_boxplot_2005_2020_by_site.png"), width = 1000, height = 600)
par(mfrow = c(1,2))
boxplot(Circumf_2005_cm ~ Site, data = growth, main = "Circumference 2005", ylab = "cm")
boxplot(Circumf_2020_cm ~ Site, data = growth, main = "Circumference 2020", ylab = "cm")
dev.off()
par(mfrow = c(1,1))

# Q9: mean 10-year growth (2010->2020) by Site
growth$Growth_10yrs <- growth$Circumf_2020_cm - growth$Circumf_2010_cm
mean10 <- aggregate(Growth_10yrs ~ Site, data = growth, FUN = mean)
sd10   <- aggregate(Growth_10yrs ~ Site, data = growth, FUN = sd)
tab10  <- merge(mean10, sd10, by = "Site", suffixes = c("_Mean", "_SD"))
write.csv(tab10, file.path(outdirT, "q09_growth10_by_site.csv"), row.names = FALSE)

# Q10: t-test of 10-year growth between sites
tt <- t.test(Growth_10yrs ~ Site, data = growth)
capture.output(tt, file = file.path(outdirT, "q10_ttest_growth10_by_site.txt"))
