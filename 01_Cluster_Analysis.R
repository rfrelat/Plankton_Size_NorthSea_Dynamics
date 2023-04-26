# 01_Cluster_Analysis.R: 
# average abundance per region and per year
# and run a quick hierarchical clustering (Fig 3)
# Last update: 25/04/2023


# 1. Load data and required functions
load("NBSS_Borner_2023.Rdata")

# Replace months by location
tab$region <- ifelse(tab$month=="9", "Buchan/Banks", "Downs")
tab$regionyear <- paste(tab$region, tab$year)

# calculate the mean abundance per year and month (=location)
# sum of abundance per year and month
sum_regionyear <- tapply(tab$orgDF, list(tab$species, tab$regionyear), sum, na.rm=TRUE)
#sumperyearmonth[is.na(sumperyearmonth)] <- 0
# number of station per year and month
n_regionyear <- tapply(tab$stationID, tab$regionyear, function(x) length(unique(x)))
# mean = sum / n
mean_regionyear <- t(sum_regionyear) / as.numeric(n_regionyear)


# Figure 3
# mean abundance is eighth rooted
heatmap(mean_regionyear**(1/8), Colv = NA, 
        col=rev(hcl.colors(25,"Spectral")))

