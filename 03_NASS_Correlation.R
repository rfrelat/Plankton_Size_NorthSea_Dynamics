

# 1. Load data and required functions ---------------
load("NBSS_Borner_2023.Rdata")
require(maps)
require(mapdata)


# 2. Set parameters --------------------------------

# Select survey
# set to 'Sep' for September survey in Buchan/Banks or 'Dec' for december in Downs 
paramM <- "Sep" #or Dec 

# Variable of interest : here abundance
varV <- "orgDF"

# Color palette
pal <- rev(hcl.colors(8,"PuOr"))

# Parameter to limit size ranges for Flowcam and Zooscan
flowS <- 14:22 
zooS <- 26:33 

# Set color per species type
# flowcam, zooscan, both
colFZ <- c("dodgerblue4", "skyblue", "slateblue2") 

# Environmental variables to be tested
keepV <- c("depth_m", "t_mean", "s_mean", "detritus", "d_shore","herr.abun")
labV <- c("depth", "temp", "sal", "turb", "dist", "her.larv") 


# 3. Data pre-processing ---------------------------

# Reduce the data to selected size classes
tab <- subset(tab, tab$typ=="flow" & tab$class%in%flowS | tab$typ=="zoo" & tab$class%in%zooS)
labZ <- paste0("F",min(flowS), max(flowS), "Z", min(zooS), max(zooS))

# Subset the data to selected region/month
if (tolower(paramM) == "sep"){
  tab <- subset(tab, month=="9")
} else {
  tab <- subset(tab, month=="12")
}


# 4. Calculate NASS ---------------------------

# Calculate the relative weights per size class
# sum abundance per size class
sumbin <- tapply(tab[,varV], tab[,"class"], sum)
# relative weight separated between Flowcam and Zooscan
relwei <- sumbin
relwei[as.character(flowS)] <- sumbin[as.character(flowS)]/sum(sumbin[as.character(flowS)])
relwei[as.character(zooS)] <- sumbin[as.character(zooS)]/sum(sumbin[as.character(zooS)])

# keep the information relative to size classes in a data.frame
infoC <- data.frame(
  "class"=as.numeric(names(relwei)),
  "relwei"=relwei,
  "mean_bins"=tab$mean_bins[match(names(relwei), tab$class)],
  "bin_range"=tab$bin_range[match(names(relwei), tab$class)],
  "n"= tapply(tab$org, tab[,"class"], sum)
)

# Compute size spectra for all data together
All_sum <- tapply(tab[,varV], tab[,"mean_bins"], sum)
# merge with information for each size class
dfull <- as.data.frame(cbind(All_sum, infoC[match(names(All_sum), infoC$mean_bins),]))
#normalization with the range of size class
dfull$All_norm <-dfull$All_sum/dfull$bin_range
#log transformation - log2 because size classes are log2 uniform
dfull$logAll <- log(dfull$All_norm, 2)
dfull$logsize <- log(dfull$mean_bins, 2)
#linear regression with relative weights
lm_full <- lm(logAll ~ logsize, weights = relwei, data = dfull)
sl_full <- as.numeric(coef(lm_full)[2])
ii_full <- as.numeric(coef(lm_full)[1])
r2_full <- summary(lm_full)$adj.r.squared


# Compute the size spectra per station
sl <- c() #slope
ii <- c() #intercept
r2 <- c() #r-square
for (i in sort(unique(tab$sampleID))){
  #select data for station i
  tabi <- subset(tab, tab$sampleID==i)
  
  # compute the sum per size class
  i_sum <- tapply(tabi[,varV], tabi[,"mean_bins"], sum)
  
  # merge with information for each size class
  di <- as.data.frame(cbind(i_sum, infoC[match(names(i_sum), infoC$mean_bins),]))
  
  #normalization with the range of size class
  di$i_norm <-di$i_sum/di$bin_range
  
  #log transformation - log2 because size classes are log2 uniform
  di$logbv <- log(di$i_norm, 2)
  di$logsize <- log(di$mean_bins, 2)
  
  #linear regression with relative weigths
  lm_size <- lm(logbv ~ logsize, weights = relwei, data = di)
  #slope
  sli <- as.numeric(coef(lm_size)[2])
  #intercept
  iii <- as.numeric(coef(lm_size)[1])
  #goodness of fit (R2)
  r2i <- summary(lm_size)$adj.r.squared
  
  sl <- c(sl, sli)
  ii <- c(ii, iii)
  r2 <- c(r2, r2i)
}

#save the NASS characteristics in a data.frame
sizespe <- data.frame(
  "sampleID"=sort(unique(tab$sampleID)),
  "slope"=sl,
  "intercept"=ii,
  "r2"=r2, 
  "N"=tapply(tab[,varV], tab$sampleID, sum)
)

# get information about the stations
sampleID <- strsplit(sizespe$sampleID, "_")
sizespe$year <- substr(sapply(sampleID, "[[", 1),7,10)
sizespe$month <- substr(sapply(sampleID, "[[", 1),5,6)
sizespe$station <- sapply(sampleID, "[[", 2)


# 5. Compare with environmental variables -----------------
# match the stations in the data.frames sizespe and env
env <- env[match(sizespe$sampleID, env$st_ID),]

# subset the environmental dataset
env1 <- env[,keepV]
names(env1) <- labV

# correlation slope vs env
comp <- data.frame(cbind("slope"=sizespe$slope, env1))
# compute correlation
corslope <- cor(comp, method = "pearson")[-1,1]

# Get the coordinates of stations
coo <- cbind(env$lon, env$lat)
# add jitter to avoid overlapping of stations
cooji <- jitter(coo, factor = 30)

# 6. Plot the figure ----------------------

# set the layout
layout(matrix(c(1,3,3,2,4,5), ncol=2), heights = c(1,1,1))

# A NASS 
par(mar=c(4,5,1,1)) #mfrow=c(3,2),
plot(logAll ~ logsize, data=dfull, 
     xlab="log2(Biovolume)", ylab="log2(Abundance)",
     col=ifelse(names(sumbin)%in%flowS,colFZ[1], colFZ[2]),
     cex=dfull$relwei*3+0.8, pch=16)

labR<-c(paste("Rsq =",round(mean(sizespe$r2),2)),
        paste("Sl =",round(mean(sizespe$slope),2)),
        paste("In =",round(mean(sizespe$intercept),2)))
legend("topright",legend=labR,bty="n")
xpred <- c(14,33)
ypred <- predict(lm_full, newdata = data.frame("logsize"=xpred))
lines(xpred, ypred)


# B histogram slope
hist(sizespe$slope, main="", xlab="slope",ylab="Frequency")
xlim <- range(sizespe$slope)
xseq <- seq(xlim[1],xlim[2], length.out = length(pal)+1)
rect(xseq[-length(xseq)],par()$usr[3],
     xseq[-1],par()$usr[3]+1,  
     col=pal, border = NA, xpd=NA)


# C Map of slopes
par(mar=c(4,4,1,1))
#all stations
plot(cooji, type="n", asp=1/cos(56*pi/180), 
     xlab="Longitude (E)", ylab="Latitude (N)") 
#add background map
map("world", fill = TRUE, 
    col="grey", add=TRUE)
# color based on PC1 score
colv <- pal[cut(sizespe$slope, xseq, include.lowest=TRUE)]
#plot the point
points(cooji,  pch=16, col=colv)


# D Time series of slopes
par(mar=c(3,4,1,1), las=1)

qsl <- tapply(sizespe$slope, sizespe$year, quantile, 
               probs=c(0.025,0.25,0.5, 0.75, 0.975))
qsl <- sapply(qsl, as.numeric)

yr <- colnames(qsl)
plot(yr, qsl[3,], ylim=range(qsl), xlab="", 
     ylab="NASS", type="n")
polygon(c(yr, rev(yr)), c(qsl[1,], rev(qsl[5,])), 
        col="grey90", border=NA)
polygon(c(yr, rev(yr)), c(qsl[2,], rev(qsl[4,])), 
        col="grey70", border=NA)
lines(yr, qsl[3,], lwd=2)

# remove plot for year 2015 (no data)
if (tolower(paramM)=="sep"){
  polygon(c(2014.2, 2015.8, 2015.8, 2014.2), c(rep(min(qsl), 2), rep(max(qsl),2)), 
          col="white", border=NA)
} else {
  polygon(c(2013.2, 2014.8, 2014.8, 2013.2), c(rep(min(qsl), 2), rep(max(qsl),2)), 
          col="white", border=NA)
}
#Add color scale
rect(par()$usr[1], xseq[-length(xseq)],
     par()$usr[1]+0.1, xseq[-1], 
     col=pal, border = NA)


# E Environmental correlation
par(mar=c(4,5,1,1), las=1)
ordE <- order(corslope)
mpc1 <- max(abs(corslope))
bke <- seq(-mpc1, mpc1, length.out=9)
cole <- pal[cut(corslope, bke, include.lowest=TRUE)]
barplot(corslope[ordE], horiz=TRUE, col=cole[ordE],
        names.arg = labV[ordE], xlab="Correlation with slope")

