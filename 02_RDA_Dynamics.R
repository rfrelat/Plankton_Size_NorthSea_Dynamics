# 02_RDADynamics.R: 
# run RDA for either Buchan/Banks in autumn 
# or Downs in winter (Fig 4 and 5)  
# Last update: 10/10/2024

# 1. Load data and required functions ---------------------
load("NBSS_Borner_2023.Rdata")
require(ade4)
require(vegan)
require(maps)
require(mapdata)
require(car)
require(nortest)

# 2. Set parameters ---------------------------------------

# Select survey
# set to 'Sep' for September survey in Buchan/Banks or 'Dec' for December in Downs 
paramM <- "Dec" #or Dec 

# Variable of interest : here abundance (orgDF), biomass (carbonDF), biovolume (biovolDF)
varV <- "orgDF"

# Color palette
pal <- hcl.colors(8,"BrBG")

# Parameter to limit size ranges for Flowcam and Zooscan
flowS <- 14:22 
zooS <- 26:33 

# Environmental variables to be tested
keepV <- c("depth_m", "t_mean", "s_mean", "detritus", "d_shore","herr.abun")
labV <- c("depth", "temp", "sal", "turb", "dist", "her.larv") 

# Set color per species type
# flowcam, zooscan, both
colFZ <- c("dodgerblue4", "skyblue", "slateblue2") 

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

# Calculate the number of individual per size class
nsize <- table(tab$class)

# Create a variable with species and size class
tab$spesize <- paste(tab$species, tab$class, sep="_")

# Compute abundance per size class, year and station
abu <- tapply(tab[,varV], list(tab$sampleID, tab$spesize), 
              sum, na.rm=TRUE)

# transform NA as 0
abu[is.na(abu)] <- 0 

dim(abu) #282 stations x 158 species and size

# keep station information
sampleID <- strsplit(rownames(abu), "_")
info <- data.frame(
  "year"=substr(sapply(sampleID, "[[", 1),7,10),
  "month"=substr(sapply(sampleID, "[[", 1),5,6),
  "station"=sapply(sampleID, "[[", 2)
)
row.names(info) <- rownames(abu)

# keep species information
colID <- strsplit(colnames(abu), "_")
infoC <- data.frame(
  "species"=sapply(colID, "[[", 1),
  "size"=sapply(colID, "[[", 2)
)
infoC$species <- as.factor(infoC$species)
infoC$size <- as.factor(infoC$size)

# match the stations in the data.frames abu and env
env <- env[match(row.names(abu), env$st_ID),]

# subset the environmental dataset
env1 <- env[,keepV]
names(env1) <- labV

# Get the coordinates
coo <- cbind(env$lon, env$lat)
# add jitter to avoid overlapping of stations
cooji <- jitter(coo, factor = 30)


# find the species type based on size (for plotting)
# 1 flowcam, 2: zooscan, 3: both 
spetyp <- table(tab$species,tab$typ)
spetypC <- apply(spetyp, 1, which.max)
spetypC[apply(spetyp>0,1,sum)==2] <- 3

clatyp <- table(tab$class,tab$typ)
clatypC <- apply(clatyp, 1, which.max)

# 4. Run the RDA --------------------------

# Reshape the 'abu' matrix to a long format for Levene's Test
abu_long <- as.data.frame(as.table(abu))

# Rename columns for clarity
colnames(abu_long) <- c("station", "spesize", "abundance")

# Convert 'station' to a factor if it's not already
abu_long$station <- as.factor(abu_long$station)

# Perform Levene's Test for homogeneity of variances
levene_test <- leveneTest(abundance ~ station, data = abu_long, center = "median")
print(levene_test)

# Perform Anderson-Darling Test for normality on abundance data
ad_test <- ad.test(abu_long$abundance)
print(ad_test)

# abundance is eighth rooted
sqabu <- abu**(1/8)

# then Hellinger transformed
abuH <- decostand(sqabu, method = "hellinger")

# Reassess homogeneity of variances after transformation
# Reshape the transformed data for Levene's Test
abuH_long <- as.data.frame(as.table(abuH))
colnames(abuH_long) <- c("station", "spesize", "abundance")
abuH_long$station <- as.factor(abuH_long$station)

# Perform Levene's Test on the transformed data
levene_test_after <- leveneTest(abundance ~ station, data = abuH_long, center = "median")
print(levene_test_after)

# # Perform Anderson-Darling Test for normality on the transformed abundance data
# ad_test_transformed <- ad.test(sqabu_long$abundance_transformed)
# print(ad_test_transformed)

# run PCA with transformed abundance
pca1 <- dudi.pca(abuH, scannf = FALSE, nf = 3)

# Detrended Correspondence Analysis
# DCA <- decorana (log1p (abuH))
# DCA

# run the RDA with PCA results and environmental variables
rda1 <- pcaiv(pca1, as.data.frame(env1),
              scannf = FALSE, nf = 3)

# ANOVA on RDA1 scores against environmental variables
p <- data.frame(RDA1 = rda1$li[,1], env1)
anova_results <- lapply(names(p)[2:ncol(p)], function(var) summary(aov(RDA1 ~ p[[var]], data = p)))
names(anova_results) <- names(p)[2:ncol(p)]
print(anova_results)

# Perform a permutation test to assess the significance of the R² value from the RDA
testR2.rda.F <- randtest(rda1, nrepet = 999)
testR2.rda.F # test of the R2 (portion of community data explained by environmental variables)
rda1$cor ## correlation with environmental variables

# labels PC1
labpc1 <- paste0("PC1: ", round(rda1$eig[1]/sum(rda1$eig)*100,1), "%",
                 " ~ ", round(rda1$eig[1]/sum(pca1$eig)*100,1), "%")

# to simplify the comparison between September and December
# we change the sign of RDA1 in December
if (paramM != "Sep"){
  #change the sign of PC1
  rda1$co[,1] <-  -rda1$co[,1]
  rda1$li[,1] <-  -rda1$li[,1]
  rda1$cor[,1] <-  -rda1$cor[,1]
}

# 5. Plot the figure ----------------------

# Define the mapping between varV values and corresponding labels
label_map <- list(orgDF = "abundance",
                  carbonDF = "biomass",
                  biovolDF = "biovolume")

# Loop through each varV value
for (v in varV) {
  # Create the plot title based on the current varV value
  title <- paste(labpc1, "- Metric:", label_map[[v]])
}

# set the layout
layout(matrix(c(1,3,3,2,4,5), ncol=2), heights = c(1.5,1,1))

# A : Species loading
par(mar=c(4,5,1,1))
# order species per average loading
ordSp <- order(tapply(rda1$co[,1], infoC$species, mean))
Spe <- factor(infoC$species, levels = levels(infoC$species)[ordSp])
colSpe <- colFZ[spetypC[match(levels(Spe),names(spetypC))]]
boxplot(rda1$co[,1]~Spe, horizontal = TRUE, col=colSpe,
        xlab="PC1", ylab="", las=1, cex.axis=0.7,
        main=title)
abline(v=0, lty=2, col="grey")
#Add color scale
xlim <- c(-1,1)*(max(abs(rda1$co[,1])))
xseq <- seq(xlim[1],xlim[2], length.out = length(pal)+1)
rect(xseq[-length(xseq)], par()$usr[3], 
     xseq[-1], par()$usr[3]+0.3,
     col=pal, border = NA)


# B : Size loading
par(mar=c(4,4,1,1))
colSiz <- colFZ[clatypC[match(levels(infoC$size),names(clatypC))]]
boxplot(rda1$co[,1]~infoC$size, horizontal = TRUE,
        col=colSiz, xlab="PC1", ylab="size class", las=1, 
        cex.axis=0.7)
rect(xseq[-length(xseq)], par()$usr[3], 
     xseq[-1], par()$usr[3]+0.3,
     col=pal, border = NA)
abline(v=0, lty=2, col="grey")


# C: Spatial distribution of the loadings
par(mar=c(4,4,1,1))
#all stations
plot(cooji, type="n", asp=1/cos(56*pi/180), 
     xlab="Longitude (E)", ylab="Latitude (N)") 
#add background map
map("world", fill = TRUE, 
    col="grey", add=TRUE)
# color based on PC1 score (absolute value scale)
mpc1 <- max(abs(rda1$li[,1]))
bkv <- seq(-mpc1, mpc1, length.out=9) 
colv <- pal[cut(rda1$li[,1], bkv, include.lowest=TRUE)]
#plot the point
points(cooji,  pch=16, col=colv)


# D: Time series loading
par(mar=c(3,4,1,1), las=1)
# calculate quantiles per year
qpc1 <- tapply(rda1$li[,1], info$year, quantile, 
               probs=c(0.025,0.25,0.5, 0.75, 0.975))
qpc1 <- sapply(qpc1, as.numeric)
yr <- colnames(qpc1)
plot(yr, qpc1[3,], ylim=range(qpc1), xlab="", 
     ylab="PC1", type="n")
polygon(c(yr, rev(yr)), c(qpc1[1,], rev(qpc1[5,])), 
        col="grey90", border=NA)
polygon(c(yr, rev(yr)), c(qpc1[2,], rev(qpc1[4,])), 
        col="grey70", border=NA)
lines(yr, qpc1[3,], lwd=2)

# remove plot for year 2015 (no data)
if (tolower(paramM)=="sep"){
  polygon(c(2014.2, 2015.8, 2015.8, 2014.2), c(rep(min(qpc1), 2), rep(max(qpc1),2)), 
          col="white", border=NA)
} else {
  polygon(c(2013.2, 2015, 2015, 2013.2), c(rep(min(qpc1), 2), rep(max(qpc1),2)), 
          col="white", border=NA)
}

#Add color scale
ylim <- c(-1,1)*(max(abs(rda1$li[,1])))
yseq <- seq(ylim[1],ylim[2], length.out = length(pal)+1)
rect(par()$usr[1], yseq[-length(yseq)],
     par()$usr[1]+0.1, yseq[-1], 
     col=pal, border = NA)


# E : Environmental loading
par(mar=c(4,5,1,1), las=1)
ordE <- order(rda1$cor[,1])
mpc1 <- max(abs(rda1$cor[,1]))
bke <- seq(-mpc1, mpc1, length.out=9)
cole <- pal[cut(rda1$cor[,1], breaks=bke, include.lowest=TRUE)]
barplot(rda1$cor[ordE,1], horiz=TRUE, col=cole[ordE],
        names.arg = labV[ordE], xlab="PC1", xlim=c(min(rda1$cor[,1])-0.2,max(rda1$cor[,1])+0.2))

rda1$cor
